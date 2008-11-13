/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "OrientationL.hpp"
#include "../2particleEventData.hpp"
#include <boost/progress.hpp>
#include "../../datatypes/vector.xml.hpp"

void
CLNOrientation::operator<<(const XMLNode& XML)
{
  XMLNode xSubNode(XML.getParentNode().getParentNode()
		   .getChildNode("ParticleData"));
  
  size_t nPart = xSubNode.nChildNode("Pt");
  
  orientationData.resize(nPart);
  
  I_cout() << "Loading orientation data....";

  boost::progress_display prog(nPart);
  
  int xml_iter = 0;

  for (size_t i = 0; i < nPart; ++i)
    {
      XMLNode xBrowseNode(xSubNode.getChildNode("Pt", &xml_iter));
      
      orientationData[i].orientation << xBrowseNode.getChildNode("U");
      orientationData[i].angularVelocity << xBrowseNode.getChildNode("O");
      ++prog;
    }
}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation";
}

void 
CLNOrientation::outputExtraPDatXML(xmlw::XmlStream& XML,
				   const CParticle& part) const
{
  XML << xmlw::tag("O")
      << orientationData[part.getID()].angularVelocity
      << xmlw::endtag("O")
      << xmlw::tag("U")
      << orientationData[part.getID()].orientation
      << xmlw::endtag("U");
}

bool 
CLNOrientation::getLineLineCollision(CPDData& PD, const Iflt& length, 
				     const CParticle& p1, const CParticle& p2) const
{ 
  // +0.1 is arbitrary to ensure a non-zero rate if angular velocities are both zero
  Iflt interpolationSize = 1.0 / (0.1 + (10.0 * orientationData[p1.getID()].angularVelocity.length()
                                              * orientationData[p2.getID()].angularVelocity.length()));

  // If interpolation size is > window size, put rate as window width
  interpolationSize = (interpolationSize > PD.dt) ? PD.dt : interpolationSize;
  
  orientationStreamType A, B;
  
  A.position = p1.getPosition();
  A.velocity = p1.getVelocity();
  A.rot.orientation = orientationData[p1.getID()].orientation;
  A.rot.angularVelocity = orientationData[p1.getID()].angularVelocity;
  
  B.position = p2.getPosition();
  B.velocity = p2.getVelocity();
  B.rot.orientation = orientationData[p2.getID()].orientation;
  B.rot.angularVelocity = orientationData[p2.getID()].angularVelocity;
  
  return recursiveRootFinder(A, B, length, interpolationSize, 0.0, PD.dt, PD.dt);
}

bool
CLNOrientation::recursiveRootFinder(orientationStreamType& A, orientationStreamType& B, const Iflt& length,
                                    const Iflt& interpolationSize, const Iflt& windowOpen, const Iflt& windowClosed, Iflt& collisionTime) const
{
  long unsigned int iter = 0;
  Iflt currentPosition = 0, x0 = 0, x1 = 0, x2 = 0, 
       upperTimeBracket, lowerTimeBracket, upperValue, 
       lowerValue, midpoint, previousMidpoint, workingValue,
       rijdotui, rijdotuj, uidotuj, alpha, beta;
  CVector<> rij, crossProduct;
  orientationStreamType rootWorkerA, rootWorkerB;
  
  while(currentPosition < windowClosed)
  {
    currentPosition = iter * interpolationSize;
    
    performRotation(A, interpolationSize);
    performRotation(B, interpolationSize);
    
    x0 = x1;
    x1 = x2;
    
    rij = A.position - B.position;
    crossProduct = A.rot.orientation.Cross(B.rot.orientation);
	
	  x2 = crossProduct % rij;
    
    // Possible root found in the x1/x2 interval
    if(std::signbit(x2) != std::signbit(x1) && iter > 0)
    {
      // Root checking routine here
      rootWorkerA = A;
      rootWorkerB = B;
      upperValue = x2;
      lowerValue = x1;
      upperTimeBracket = currentPosition;
      lowerTimeBracket = currentPosition - interpolationSize;
      previousMidpoint = currentPosition;
      midpoint = currentPosition;
      
      while(std::fabs(upperTimeBracket - lowerTimeBracket) > eps)
      {
        // Bisection root finding
        previousMidpoint = midpoint;
        midpoint = (lowerTimeBracket + upperTimeBracket) / 2;
        
        performRotation(rootWorkerA, midpoint - currentPosition);
        performRotation(rootWorkerB, midpoint - currentPosition);
        
        rij = rootWorkerA.position - rootWorkerB.position;
        crossProduct = rootWorkerA.rot.orientation.Cross(rootWorkerB.rot.orientation);
        
        workingValue = crossProduct % rij;
        
        if(std::signbit(workingValue) == std::signbit(lowerValue))
        {
          lowerTimeBracket = midpoint;
        }
        else
        {
          upperTimeBracket = midpoint;
        }
      }
      
      // Now upperTimeBracket will give us the collision time.
      // Check that the root is valid.
      rijdotui = rij % rootWorkerA.rot.orientation;
      rijdotuj = rij % rootWorkerB.rot.orientation;
      uidotuj = rootWorkerA.rot.orientation % rootWorkerB.rot.orientation;
      
      alpha = -1.0 * (rijdotui - (rijdotuj * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
      beta  = -1.0 * (rijdotuj - (rijdotui * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
      
      if(std::fabs(alpha) < (length/2) && std::fabs(beta) < (length/2))
      {
        collisionTime = upperTimeBracket;
        return true;
      }
    }
    
    // Conditions for interpolating:
    // - Sign is different to what extrapoling the previous segment indicated
    // - The new point is still the same sign as the old one
    if(iter > 1 && (std::signbit(x2) != std::signbit(x1 + x1 - x0)) && std::signbit(x2) == std::signbit(x1))
    {
      performRotation(A, -1.0 * interpolationSize);
      performRotation(B, -1.0 * interpolationSize);
      
      if(recursiveRootFinder(A, B, length, interpolationSize/10.0, currentPosition - interpolationSize, currentPosition, collisionTime))
      {
        return true;
      }
    }
    
    iter++;
  }
  
  return false;
}

C2ParticleData 
CLNOrientation::runLineLineCollision(const CIntEvent&) const
{ D_throw() << "Not implemented"; }

void
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  orientationStreamType ostPart;
  
  ostPart.velocity = part.getVelocity();
  ostPart.position = part.getPosition();
  ostPart.rot.angularVelocity = orientationData[part.getID()].angularVelocity;
  ostPart.rot.orientation = orientationData[part.getID()].orientation;
  
  performRotation(ostPart, dt);
  
  part.getPosition() =  ostPart.position;
  orientationData[part.getID()].orientation = ostPart.rot.orientation;
}

void
CLNOrientation::performRotation(orientationStreamType& osret, const Iflt& dt) const
{
  if(NDIM != 3) { D_throw() << "Implemented only for 3D rotations"; }
  else
  {
    // Linear dynamics
    osret.position += (osret.velocity * dt);
    
    Iflt angle = osret.rot.angularVelocity.length() * dt;
  
    Iflt v1 = osret.rot.angularVelocity.unitVector()[0];
    Iflt v2 = osret.rot.angularVelocity.unitVector()[1];
    Iflt v3 = osret.rot.angularVelocity.unitVector()[2];
  
    // axis is not undefined and angle is not zero
    if(!(v1 == 0 && v2 == 0 && v3 == 0) && angle != 0)
    {	
    
      Iflt matrix[3][3];
    
      Iflt cos_term = cos(angle);
      Iflt sin_term = sin(angle);
  
      matrix[0][0] = pow(v1, 2) + (pow(v2 ,2) + pow(v3, 2))*(cos_term);
      matrix[0][1] = (v1 * v2 * (1 - cos_term)) - (v3 * sin_term);
      matrix[0][2] = (v1 * v3 * (1 - cos_term)) + (v2 * sin_term);
      matrix[1][0] = (v1 * v2 * (1 - cos_term)) + (v3 * sin_term);
      matrix[1][1] = pow(v2, 2) + (pow(v3 ,2) + pow(v1, 2))*(cos_term);
      matrix[1][2] = (v2 * v3 * (1 - cos_term)) - (v1 * sin_term);
      matrix[2][0] = (v3 * v1 * (1 - cos_term)) - (v2 * sin_term);
      matrix[2][1] = (v2 * v3 * (1 - cos_term)) + (v1 * sin_term);
      matrix[2][2] = pow(v3, 2) + (pow(v1 ,2) + pow(v2, 2))*(cos_term);
    
      CVector<> tempvec;
      tempvec[0] = (matrix[0][0] * osret.rot.orientation[0]) + (matrix[0][1] * osret.rot.orientation[1]) + (matrix[0][2] * osret.rot.orientation[2]);
      tempvec[1] = (matrix[1][0] * osret.rot.orientation[0]) + (matrix[1][1] * osret.rot.orientation[1]) + (matrix[1][2] * osret.rot.orientation[2]);
      tempvec[2] = (matrix[2][0] * osret.rot.orientation[0]) + (matrix[2][1] * osret.rot.orientation[1]) + (matrix[2][2] * osret.rot.orientation[2]);
    
      osret.rot.orientation = tempvec;
    }
  }
}

C1ParticleData 
CLNOrientation::runAndersenWallCollision(const CParticle& part, 
			 const CVector<>& vNorm,
			 const Iflt& sqrtT
			 ) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";
}
  
C1ParticleData 
CLNOrientation::randomGaussianEvent(const CParticle& part, 
				    const Iflt& sqrtT) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";  
}

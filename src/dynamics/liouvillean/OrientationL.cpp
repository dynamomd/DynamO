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
CLNOrientation::getLineLineCollision(const CPDData& PD, const Iflt& length, 
				     const CParticle& p1, const CParticle& p2,
				     const Iflt& twindow
				     ) const
{ 
  // +0.1 is arbitrary to ensure a non-zero rate if angular velocities are both zero
  Iflt interpolationSize = 1.0 / (0.1 + (10.0 * orientationData[p1.getID()].angularVelocity.length()
                                              * orientationData[p2.getID()].angularVelocity.length()));

  // If interpolation size is > window size, put rate as window width
  interpolationSize = (interpolationSize > twindow) ? twindow : interpolationSize;
  
  return recursiveRootFinder(interpolationSize, 0.0, twindow);
}

bool
CLNOrientation::recursiveRootFinder(const Iflt& interpolationSize, const Iflt& window_open, const Iflt& window_closed) const
{ D_throw() << "Not implemented"; }

C2ParticleData 
CLNOrientation::runLineLineCollision(const CIntEvent&) const
{ D_throw() << "Not implemented"; }

void
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  orientationStreamReturnType osret = performRotation(part, dt);
  
  CLNewton::streamParticle(part, dt);
  
  orientationData[part.getID()].orientation = osret.rot.orientation;
}

CLNOrientation::orientationStreamReturnType
CLNOrientation::performRotation(CParticle& part, const Iflt& dt) const
{
  if(NDIM != 3) { D_throw() << "Implemented only for 3D rotations"; }
  else
  {
    orientationStreamReturnType osret;
    
    // Linear dynamics
    osret.velocity = part.getVelocity();
    osret.position = part.getPosition() + (osret.velocity * dt);
    osret.rot.angularVelocity = orientationData[part.getID()].angularVelocity;
    osret.rot.orientation = orientationData[part.getID()].orientation;
  
    // Angular dynamics
    osret.rot.angularVelocity = orientationData[part.getID()].angularVelocity;
    
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
  
    return osret;
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

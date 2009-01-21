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
  XMLNode xSubNode = XML.getChildNode("OrientationData");
  
  size_t nPart = xSubNode.nChildNode("P");
    
  orientationData.resize(nPart);
  
  I_cout() << "Loading orientation data....";

  boost::progress_display prog(nPart);
  
  int xml_iter = 0;

  for (size_t i = 0; i < nPart; ++i)
    {
      XMLNode xBrowseNode(xSubNode.getChildNode("P", &xml_iter));
      
      orientationData[i].orientation << xBrowseNode.getChildNode("U");
      orientationData[i].angularVelocity << xBrowseNode.getChildNode("O");
      
      Iflt oL = orientationData[i].orientation.length();

      if (!(oL > 0.0))
	D_throw() << "Particle ID " << i 
		  << " orientation vector is zero!";

      //Makes the vector a unit vector
      orientationData[i].orientation /= oL;
      
      ++prog;
    }
}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation"
      << xmlw::tag("OrientationData");
  
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    XML << xmlw::tag("P") 
	<< xmlw::attr("ID") << part.getID()
	<< xmlw::tag("O") 
	<< orientationData[part.getID()].angularVelocity
	<< xmlw::endtag("O")
	<< xmlw::tag("U")
	<< orientationData[part.getID()].orientation
	<< xmlw::endtag("U")
	<< xmlw::endtag("P");

  XML << xmlw::endtag("OrientationData");
}

bool 
CLNOrientation::getLineLineCollision(CPDData& PD, const Iflt& length, 
				     const CParticle& p1, const CParticle& p2) const
{ 
  Iflt windowSize = PD.dt;  
  Iflt currentClock = 0.0;

  Iflt t_low = 0.0;
  Iflt t_up = PD.dt;

  Iflt t_low_working = t_low;
  Iflt t_up_working = t_up;

  Iflt f0, f1, f2, root1, root2, foundRoot, initialRoot;
  
  // Set up pair of lines as passable objects
  orientationStreamType A, B;

  A.position = PD.rij;
  A.velocity = PD.vij;
  A.rot.orientation = orientationData[p1.getID()].orientation;
  A.rot.angularVelocity = orientationData[p1.getID()].angularVelocity;
  
  B.position = CVector<>(0);
  B.velocity = CVector<>(0);
  B.rot.orientation = orientationData[p2.getID()].orientation;
  B.rot.angularVelocity = orientationData[p2.getID()].angularVelocity;

  // Get Frenkel second derivative maximum
  Iflt maxSecondDeriv = F_secondDeriv_max(A, B, length);

  while(t_low < t_up)
  {
    // Calculate f0, f1, f2 at t = t_low
    f0 = F_zeroDeriv(A, B);
    f1 = F_firstDeriv(A, B);
    f2 = F_secondDeriv(A, B);

    // Get a root
    initialRoot = quadraticSolution(ROOT_SMALLEST_POSITIVE, f0, f1, f2);

    if(initialRoot > t_up)
    {
      D_throw() << "Root guess is outside boundary";
      continue;
    }

    // Reset the lower search limit
    t_low_working += quadraticSolution(ROOT_SMALLEST_POSITIVE, f0, f1, maxSecondDeriv * ((f0 < 0) ? 0.5 : -0.5));
    
    // FoundRoot now has the new root
    foundRoot = quadraticRootFinder(A, B, root1);

    if(foundRoot > t_up)
    {
      D_throw() << "Root found is greater than t_up";
      continue;
    }    

    // wind the lines to the new root position, and shrink the t_up value
    performRotation(A, foundRoot);
    performRotation(B, foundRoot);

    t_up_working = foundRoot - ((2 * fabs(F_firstDeriv(A, B)))/maxSecondDeriv);
    
    // wind the lines to the t_low position and search again
    performRotation(A, t_low_working - foundRoot);
  }
  
  return false;
  //return recursiveRootFinder(A, B, length, interpolationSize, 0.0, PD.dt, PD.dt);
}

bool
CLNOrientation::recursiveRootFinder(orientationStreamType& A, orientationStreamType& B, const Iflt& length,
                                    const Iflt& interpolationSize, const Iflt& windowOpen, const Iflt& windowClosed, Iflt& collisionTime) const
{
  /*
  long unsigned int iter = 0;
  Iflt currentPosition = windowOpen, x0 = 0, x1 = 0, x2 = 0, 
       upperTimeBracket, lowerTimeBracket, upperValue, 
       lowerValue, midpoint, previousMidpoint, workingValue;
  CVector<> rij, crossProduct;
  orientationStreamType rootWorkerA, rootWorkerB;
  
  while(currentPosition < windowClosed)
  {
    x0 = x1;
    x1 = x2;
    
    rij = A.position - B.position;
    crossProduct = A.rot.orientation.Cross(B.rot.orientation);
	
	  x2 = crossProduct % rij;
    
    // Possible root found in the x1/x2 interval
    if(std::signbit(x2) != std::signbit(x1) && iter > 1)
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
      
      while(std::fabs(upperTimeBracket - lowerTimeBracket) >  (eps*lowerTimeBracket))
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
      collisionPoints cp = getCollisionPoints(rootWorkerA, rootWorkerB);
      
      if(std::fabs(cp.alpha) < (length/2) && std::fabs(cp.beta) < (length/2))
      {
        collisionTime = (upperTimeBracket + lowerTimeBracket) / 2.0;
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
    performRotation(A, interpolationSize);
    performRotation(B, interpolationSize);
    currentPosition += interpolationSize;
  }
  */
  return false;
}

Iflt
CLNOrientation::quadraticSolution(const int returnType, Iflt A, Iflt B, Iflt C) const
{
  Iflt discriminant = (B * B) - (4 * A * C);

  if(discriminant < 0)
  {
    D_throw() << "Determinant of less than zero returned";
    return 0;
  }
  
  Iflt root1 = ((-1.0 * B) + sqrt(discriminant)) / (2 * A);
  Iflt root2 = ((-1.0 * B) - sqrt(discriminant)) / (2 * A);

  if(returnType == ROOT_SMALLEST_EITHER)
  {
    return (fabs(root1) < fabs(root2)) ? root1 : root2;
  }

  else if(returnType == ROOT_LARGEST_EITHER)
  {
    return (fabs(root1) < fabs(root2)) ? root2 : root1;
  }
  else
  {    
    if(root1 > 0 && root2 > 0) // Both roots positive
    {
      switch(returnType)
      {
        case ROOT_LARGEST_NEGATIVE:
        case ROOT_SMALLEST_NEGATIVE:
          D_throw() << "Both roots positive";
          return 0;
          break;
        case ROOT_SMALLEST_POSITIVE:
          return ((root1 < root2) ? root1 : root2);
          break;
        case ROOT_LARGEST_POSITIVE:
          return ((root1 > root2) ? root1 : root2);
	  break;
      }
    }
    else if(root1 < 0 && root2 < 0) // Both roots negative
    {
      switch(returnType)
      {
        case ROOT_LARGEST_POSITIVE:
        case ROOT_SMALLEST_POSITIVE:
          D_throw() << "Both roots negative";
          return 0;
          break;
        case ROOT_SMALLEST_NEGATIVE:
          return ((root1 > root2) ? root1 : root2);
          break;
        case ROOT_LARGEST_NEGATIVE:
          return ((root1 < root2) ? root1 : root2);
          break;
      }
    }
    else // Roots are different signs
    {
      switch(returnType)
      {
        case ROOT_LARGEST_POSITIVE:
        case ROOT_SMALLEST_POSITIVE:
          return ((root1 > root2) ? root1 : root2);
          break;
        case ROOT_LARGEST_NEGATIVE:
        case ROOT_SMALLEST_NEGATIVE:
          return ((root1 < root2) ? root1 : root2);
          break;
      }
    }
  } 
}

Iflt
CLNOrientation::F_zeroDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;
  return ((A.rot.orientation.Cross(B.rot.orientation)) % deltaR);
}

Iflt 
CLNOrientation::F_firstDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;
  CVector<> deltaW = A.rot.angularVelocity - B.rot.angularVelocity;
  CVector<> deltaV = A.velocity - B.velocity;
  
  return (
    (A.rot.orientation % deltaR) * (deltaW % B.rot.orientation)
  ) + (
    (B.rot.orientation % deltaR) * (deltaW % A.rot.orientation)
  ) - (
    (deltaW % deltaR) * (A.rot.orientation % B.rot.orientation)
  ) + ( 
    (A.rot.orientation.Cross(B.rot.orientation) % deltaV)
  );
}

Iflt
CLNOrientation::F_firstDeriv_max(orientationStreamType A, orientationStreamType B, const Iflt length) const
{
  Iflt absDeltaW = (A.rot.angularVelocity - B.rot.angularVelocity).length();
  Iflt absDeltaV = (A.velocity - B.velocity).length();

  return ((length * absDeltaW) + absDeltaV);
}

Iflt 
CLNOrientation::F_secondDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;
  CVector<> deltaW = A.rot.angularVelocity - B.rot.angularVelocity;
  CVector<> deltaV = A.velocity - B.velocity;
  
  return 2.0 * (
    (
      (A.rot.orientation % deltaV) * (deltaW % B.rot.orientation)
    ) + (
      (B.rot.orientation % deltaV) * (deltaW % A.rot.orientation)
    ) - (
      (A.rot.orientation % B.rot.orientation) * (deltaW % deltaV)
    )
  ) - (
    (deltaW % deltaR) * (deltaW % (A.rot.orientation.Cross(B.rot.orientation)))
  ) + (
    (A.rot.orientation % deltaR) * (B.rot.orientation % (A.rot.angularVelocity.Cross(B.rot.angularVelocity)))
  ) + (
    (B.rot.orientation % deltaR) * (A.rot.orientation % (A.rot.angularVelocity.Cross(B.rot.angularVelocity)))
  ) + (
    (deltaW % A.rot.orientation) * (deltaR % (B.rot.angularVelocity.Cross(B.rot.orientation)))
  ) + (
    (deltaW % B.rot.orientation) * (deltaR % (A.rot.angularVelocity.Cross(A.rot.orientation)))
  ); 
}

Iflt
CLNOrientation::F_secondDeriv_max(orientationStreamType A, orientationStreamType B, Iflt length) const
{
  Iflt absDeltaW = (A.rot.angularVelocity - B.rot.angularVelocity).length();
  Iflt absDeltaV = (A.velocity - B.velocity).length();

  return absDeltaW * ((2 * absDeltaV) + (length * (A.rot.angularVelocity.length() + B.rot.angularVelocity.length())));
}


C2ParticleData 
CLNOrientation::runLineLineCollision(const CIntEvent& eevent, const Iflt& length) const
{
  updateParticlePair(eevent.getParticle1(), eevent.getParticle2());

  C2ParticleData retVal(eevent.getParticle1(), eevent.getParticle2(),
                        Sim->Dynamics.getSpecies(eevent.getParticle1()),
                        Sim->Dynamics.getSpecies(eevent.getParticle2()),
                        CORE);
  
  Sim->Dynamics.BCs().setPBC(retVal.rij, retVal.vijold);
  // End copied section
  
  CVector<> vr, uPerp, u1dot, u2dot;
  orientationStreamType A, B;
  
  // Assume lines are the same mass... use the mass of line #1
  Iflt mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt inertia = (mass * length * length)/12.0;
  
  A.position = retVal.particle1_.getParticle().getPosition();
  A.velocity = retVal.particle1_.getParticle().getVelocity();
  A.rot.angularVelocity = orientationData[retVal.particle1_.getParticle().getID()].angularVelocity;
  A.rot.orientation = orientationData[retVal.particle1_.getParticle().getID()].orientation;
  
  B.position = retVal.particle2_.getParticle().getPosition();
  B.velocity = retVal.particle2_.getParticle().getVelocity();
  B.rot.angularVelocity = orientationData[retVal.particle2_.getParticle().getID()].angularVelocity;
  B.rot.orientation = orientationData[retVal.particle2_.getParticle().getID()].orientation;
  
  uPerp = A.rot.orientation.Cross(B.rot.orientation).unitVector();
  
  collisionPoints cp = getCollisionPoints(A, B);
  
  u1dot = A.rot.angularVelocity.Cross(A.rot.orientation) * cp.alpha;
  u2dot = B.rot.angularVelocity.Cross(B.rot.orientation) * cp.beta;
  
  vr = (A.velocity - B.velocity) + u1dot - u2dot;
  
  Iflt alpha = -1.0 * (vr % uPerp);
  alpha /= ((1.0/mass) + ((std::pow(cp.alpha, 2) + std::pow(cp.beta, 2))/(2.0 * inertia)));
  
  retVal.rvdot = retVal.rij % retVal.vijold;
  retVal.dP = uPerp * alpha;
  
  const_cast<CParticle&>(eevent.getParticle1()).getVelocity() += retVal.dP / mass;
  const_cast<CParticle&>(eevent.getParticle2()).getVelocity() -= retVal.dP / mass;
  
  orientationData[eevent.getParticle1().getID()].angularVelocity = A.rot.angularVelocity - ((A.rot.orientation* (cp.alpha / inertia)).Cross(retVal.dP));
  orientationData[eevent.getParticle2().getID()].angularVelocity = B.rot.angularVelocity + ((B.rot.orientation* (cp.beta / inertia)).Cross(retVal.dP));

  return retVal;
}

CLNOrientation::collisionPoints
CLNOrientation::getCollisionPoints(orientationStreamType& A, orientationStreamType& B) const
{
  collisionPoints retVal;
  CVector<> rij;
  Iflt rijdotui, rijdotuj, uidotuj;
  
  rij = A.position - B.position;
  rijdotui = rij % A.rot.orientation;
  rijdotuj = rij % B.rot.orientation;
  uidotuj = A.rot.orientation % B.rot.orientation;
    
  retVal.alpha = -1.0 * (rijdotui - (rijdotuj * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
  retVal.beta  = -1.0 * (rijdotuj - (rijdotui * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
  
  return retVal;
}

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

Iflt
CLNOrientation::quadraticRootFinder(orientationStreamType A, orientationStreamType B, Iflt initialJump) const
{
  Iflt currentValue = 1.0, currentTime = initialJump;
  Iflt root = 0.0;
  Iflt f0 = 0.0, f1 = 0.0, f2 = 0.0;

  performRotation(A, initialJump);
  performRotation(B, initialJump);

  while(currentValue != 0.0)
  {
    f0 = F_zeroDeriv(A, B);
    f1 = F_firstDeriv(A, B);
    f2 = F_secondDeriv(A, B);

    currentValue = f0;

    root = quadraticSolution(ROOT_SMALLEST_EITHER, f0, f1, f2);
   
    currentTime += root;
    performRotation(A, root);
    performRotation(B, root);
  }

  return currentTime;
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

void 
CLNOrientation::initLineOrientations(const Iflt& length)
{
  orientationData.resize(Sim->vParticleList.size());
  
  I_cout() << "Initialising the line orientations";

  Iflt factor = std::sqrt(12.0/(length*length));

  CVector<> angVelCrossing;

  for (size_t i = 0; i < Sim->vParticleList.size(); ++i)
    {      
      //Assign the new velocities
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        orientationData[i].orientation[iDim] = Sim->normal_sampler();

      orientationData[i].orientation = orientationData[i].orientation.unitVector();

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
        angVelCrossing[iDim] = Sim->normal_sampler();
      }
      
      orientationData[i].angularVelocity 
	= orientationData[i].orientation.Cross(angVelCrossing).unitVector() 
	* Sim->normal_sampler() * factor;
    }  
}

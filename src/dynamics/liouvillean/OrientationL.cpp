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
#include "../units/units.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filter/base64.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/stream_sink.hpp>
#include <boost/iostreams/device/stream_source.hpp>
#include <boost/iostreams/filter/base64cleaner.hpp>
#include <boost/iostreams/filter/linewrapout.hpp>
#include "../../extcode/mathtemplates.hpp"


void 
CLNOrientation::initialise() 
{
  CLiouvillean::initialise();

  Iflt sumEnergy(0.0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      sumEnergy += Sim->Dynamics.getSpecies(part).getMass()
	* orientationData[part.getID()].angularVelocity.square();
    }
  
  sumEnergy *= 0.5 * Sim->Dynamics.units().unitLength() 
    * Sim->Dynamics.units().unitLength()
    / (12.0 * Sim->Dynamics.units().unitEnergy());

  I_cout() << "System Rotational Energy " << sumEnergy
	   << "\nRotational kT " << sumEnergy / Sim->lN;
}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation";
}

bool 
CLNOrientation::getLineLineCollision(CPDData& PD, const Iflt& length, 
				     const CParticle& p1, const CParticle& p2) const
{  
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(p1))
    D_throw() << "Particle1 " << p1.getID() << " is not up to date";

  if (!isUpToDate(p2))
    D_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

  Iflt t_low = 0.0;
  Iflt t_high = PD.dt;
  
  // Set up pair of lines as passable objects
  orientationStreamType A(PD.rij, PD.vij, orientationData[p1.getID()].orientation, 
			  orientationData[p1.getID()].angularVelocity),
    B(CVector<>(0), CVector<>(0), 
      orientationData[p2.getID()].orientation, 
      orientationData[p2.getID()].angularVelocity);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    //Shift the lower bound up so we don't find the same root again
    t_low += fabs(2.0 * F_firstDeriv(A, B))
      / F_secondDeriv_max(A, B, length);
  
  //Find window delimited by discs
  IfltPair dtw = discIntersectionWindow(A, B, length);
  
  if(dtw.alpha > t_low)
    t_low = dtw.alpha;
  
  if(dtw.beta < t_high)
    t_high = dtw.beta;
  
  Iflt root = frenkelRootSearch(A, B, length, t_low, t_high);

  if (root != HUGE_VAL) 
    { 
      PD.dt = root; 
      return true; 
    }
  else 
    return false;
}


Iflt
CLNOrientation::frenkelRootSearch(const orientationStreamType A, 
				  const orientationStreamType B, 
				  Iflt length, Iflt t_low, Iflt t_high) const
{  
  Iflt root = 0.0;
	
  while(t_high > t_low)
    {
      root = quadraticRootHunter(A, B, length, t_low, t_high);

      if (root == HUGE_VAL) return HUGE_VAL;
      
      Iflt temp_high = t_high;
      do {
	// Artificial boundary just below root
	orientationStreamType tempA(A), tempB(B);
	
	performRotation(tempA, root);
	performRotation(tempB, root);
	
	Iflt Fdoubleprimemax = F_secondDeriv_max(tempA, tempB, length);
	
	temp_high = root - (fabs(2.0 * F_firstDeriv(tempA, tempB))
			    / Fdoubleprimemax);
	
	if ((temp_high < t_low) || (Fdoubleprimemax == 0)) break;
	
	Iflt temp_root = quadraticRootHunter(A, B, length, t_low, temp_high);
	
	if (temp_root == HUGE_VAL) 
	  break;
	else 
	    root = temp_root;
	
      } while(temp_high > t_low);
      
      // At this point $root contains earliest valid root guess.
      // Check root validity.
      orientationStreamType tempA(A), tempB(B);
      performRotation(tempA, root);
      performRotation(tempB, root);
      
      IfltPair cp = getCollisionPoints(tempA, tempB);
      
      if(fabs(cp.alpha) < length/2.0 && fabs(cp.beta) < length/2.0)
        return root;
      else
        t_low = root + ((2.0 * fabs(F_firstDeriv(tempA, tempB)))
			/ F_secondDeriv_max(tempA, tempB, length));
    }
    
  return HUGE_VAL;
}


bool
CLNOrientation::quadraticSolution(Iflt& returnVal, const int returnType, 
				  Iflt C, Iflt B, Iflt A) const
{
  Iflt root1(0), root2(0);
  // Contingency: if A = 0, not a quadratic = linear
  if(A == 0)
    {
      if(B == 0) return false;
      
      root1 = -1.0 * C / B;
      root2 = root1;
    }
  else
    {
      Iflt discriminant = (B * B) - (4 * A * C);
    
      if (discriminant < 0) return false;
    
      //This avoids a cancellation of errors. See
      //http://en.wikipedia.org/wiki/Quadratic_equation#Floating_point_implementation
      Iflt t((B < 0)
	     ? -0.5 * (B-sqrt(discriminant))
	     : -0.5 * (B+sqrt(discriminant)));
    
      root1 = t / A;
      root2 = C / t;
    }

  if(returnType == ROOT_SMALLEST_EITHER)
    {
      returnVal = (fabs(root1) < fabs(root2)) ? root1 : root2;
      return true;
    }
  else if(returnType == ROOT_LARGEST_EITHER)
    {
      returnVal = (fabs(root1) < fabs(root2)) ? root2 : root1;
      return true;
    }
  else
    {    
      if(root1 > 0 && root2 > 0) // Both roots positive
	{
	  switch(returnType)
	    {
	    case ROOT_LARGEST_NEGATIVE:
	    case ROOT_SMALLEST_NEGATIVE:
	      //I_cerr() << "Both roots positive";
	      return false;
	      break;
	    case ROOT_SMALLEST_POSITIVE:
	      returnVal = ((root1 < root2) ? root1 : root2);
	      return true;
	      break;
	    case ROOT_LARGEST_POSITIVE:
	      returnVal = ((root1 > root2) ? root1 : root2);
	      return true;
	      break;
	    }
	}
      else if(root1 < 0 && root2 < 0) // Both roots negative
	{
	  switch(returnType)
	    {
	    case ROOT_LARGEST_POSITIVE:
	    case ROOT_SMALLEST_POSITIVE:
	      //I_cerr() << "Both roots negative";
	      return false;
	      break;
	    case ROOT_SMALLEST_NEGATIVE:
	      returnVal = ((root1 > root2) ? root1 : root2);
	      return true;
	      break;
	    case ROOT_LARGEST_NEGATIVE:
	      returnVal = ((root1 < root2) ? root1 : root2);
	      return true;
	      break;
	    }
	}
      else // Roots are different signs
	{
	  switch(returnType)
	    {
	    case ROOT_LARGEST_POSITIVE:
	    case ROOT_SMALLEST_POSITIVE:
	      returnVal = ((root1 > root2) ? root1 : root2);
	      return true;
	      break;
	    case ROOT_LARGEST_NEGATIVE:
	    case ROOT_SMALLEST_NEGATIVE:
	      returnVal = ((root1 < root2) ? root1 : root2);
	      return true;
	      break;
	    }
	}
    }
  
  D_throw() << "Unexpected end-of-function reached.  Did you specify a valid root type?";
  return false;
}

Iflt
CLNOrientation::F_zeroDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;
  return ((A.orientation.Cross(B.orientation)) % deltaR);
}

Iflt 
CLNOrientation::F_firstDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;
  CVector<> deltaW = A.angularVelocity - B.angularVelocity;
  CVector<> deltaV = A.velocity - B.velocity;
  
  return 
    ((A.orientation % deltaR) * (deltaW % B.orientation)) + 
    ((B.orientation % deltaR) * (deltaW % A.orientation)) - 
    ((deltaW % deltaR) * (A.orientation % B.orientation)) + 
    ((A.orientation.Cross(B.orientation) % deltaV));
}

Iflt
CLNOrientation::F_firstDeriv_max(orientationStreamType A, orientationStreamType B, 
				 const Iflt length) const
{
  Iflt absDeltaW = (A.angularVelocity - B.angularVelocity).length();
  Iflt absDeltaV = (A.velocity - B.velocity).length();

  return ((length * absDeltaW) + absDeltaV);
}

Iflt 
CLNOrientation::F_secondDeriv(orientationStreamType A, orientationStreamType B) const
{
  CVector<> deltaR = A.position - B.position;  
  CVector<> deltaW = A.angularVelocity - B.angularVelocity;
  CVector<> deltaV = A.velocity - B.velocity;
  
  return 
    2.0 * (
	   ((A.orientation % deltaV) * (deltaW % B.orientation)) + 
	   ((B.orientation % deltaV) * (deltaW % A.orientation)) - 
	   ((A.orientation % B.orientation) * (deltaW % deltaV))
	   ) - 
    ((deltaW % deltaR) * (deltaW % (A.orientation.Cross(B.orientation)))) + 
    ((A.orientation % deltaR) * (B.orientation % (A.angularVelocity.Cross(B.angularVelocity)))) + 
    ((B.orientation % deltaR) * (A.orientation % (A.angularVelocity.Cross(B.angularVelocity)))) + 
    ((deltaW % A.orientation) * (deltaR % (B.angularVelocity.Cross(B.orientation)))) + 
    ((deltaW % B.orientation) * (deltaR % (A.angularVelocity.Cross(A.orientation)))); 
}

Iflt
CLNOrientation::F_secondDeriv_max(orientationStreamType A, orientationStreamType B, Iflt length) const
{
  Iflt absDeltaW = (A.angularVelocity - B.angularVelocity).length();
  Iflt absDeltaV = (A.velocity - B.velocity).length();

  return absDeltaW * ((2 * absDeltaV) + (length * (A.angularVelocity.length() + B.angularVelocity.length())));
}


C2ParticleData 
CLNOrientation::runLineLineCollision(const CIntEvent& eevent, const Iflt& elasticity, const Iflt& length) const
{
  const CParticle& particle1 = Sim->vParticleList[eevent.getParticle1ID()];
  const CParticle& particle2 = Sim->vParticleList[eevent.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  C2ParticleData retVal(particle1, particle2,
                        Sim->Dynamics.getSpecies(particle1),
                        Sim->Dynamics.getSpecies(particle2),
                        CORE);
  
  Sim->Dynamics.BCs().setPBC(retVal.rij, retVal.vijold);

  retVal.rvdot = retVal.rij % retVal.vijold;
  
  orientationStreamType A(retVal.rij, retVal.vijold, orientationData[particle1.getID()].orientation, 
			  orientationData[particle1.getID()].angularVelocity),
    B(CVector<>(0), CVector<>(0), orientationData[particle2.getID()].orientation, 
      orientationData[particle2.getID()].angularVelocity);
  
  CVector<> uPerp = A.orientation.Cross(B.orientation).unitVector();
  
  IfltPair cp = getCollisionPoints(A, B);

  // \Delta {\bf v}_{imp}
  CVector<> vr = (A.velocity - B.velocity) 
    + (cp.alpha * A.angularVelocity.Cross(A.orientation)) 
    - (cp.beta * B.angularVelocity.Cross(B.orientation));
  
  Iflt mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt inertia = (mass * length * length) / 12.0;

  retVal.dP = uPerp 
    * (((vr % uPerp) * (1.0 + elasticity)) 
       / ((2.0/mass) + ((cp.alpha * cp.alpha + cp.beta * cp.beta)/inertia)));  
  
  const_cast<CParticle&>(particle1).getVelocity() -= retVal.dP / mass;
  const_cast<CParticle&>(particle2).getVelocity() += retVal.dP / mass;

  orientationData[particle1.getID()].angularVelocity 
    -= (cp.alpha / inertia) * (A.orientation.Cross(retVal.dP));

  orientationData[particle2.getID()].angularVelocity 
    += (cp.beta / inertia) * (B.orientation.Cross(retVal.dP));

  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->dSysTime;

  return retVal;
}

CLNOrientation::IfltPair
CLNOrientation::getCollisionPoints(orientationStreamType& A, orientationStreamType& B) const
{
  IfltPair retVal;
  
  CVector<> rij = A.position - B.position;
  Iflt rijdotui = rij % A.orientation;
  Iflt rijdotuj = rij % B.orientation;
  Iflt uidotuj = A.orientation % B.orientation;
    
  retVal.alpha = - (rijdotui - (rijdotuj * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
  retVal.beta  = (rijdotuj - (rijdotui * uidotuj)) / (1.0 - std::pow(uidotuj, 2));
  
  return retVal;
}

void
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  orientationStreamType ostPart(part.getPosition(), part.getVelocity(), 
				orientationData[part.getID()].orientation, 
				orientationData[part.getID()].angularVelocity);

  performRotation(ostPart, dt);
  
  part.getPosition() =  ostPart.position;
  orientationData[part.getID()].orientation = ostPart.orientation;
}

void
CLNOrientation::performRotation(orientationStreamType& osret, const Iflt& dt) const
{

  if (NDIM != 3) { D_throw() << "Implemented only for 3D rotations"; }

  else
    {
      // Linear dynamics
      osret.position += (osret.velocity * dt);
    
      Iflt angle = osret.angularVelocity.length();
  
      CVector<> v = osret.angularVelocity / angle;
      angle *= dt;
      CVector<> vsq(v);
      
      for (size_t i = 0; i < NDIM; ++i)
	vsq[i] *= vsq[i];

      // axis is not undefined and angle is not zero
      if(!(v[0] == 0 && v[1] == 0 && v[2] == 0) && angle != 0)
	{	
    
	  Iflt matrix[3][3];
    
	  Iflt cos_term = cos(angle);
	  Iflt sin_term = sin(angle);
  
	  matrix[0][0] = vsq[0] + (vsq[1] + vsq[2])*(cos_term);
	  matrix[0][1] = (v[0] * v[1] * (1 - cos_term)) - (v[2] * sin_term);
	  matrix[0][2] = (v[0] * v[2] * (1 - cos_term)) + (v[1] * sin_term);
	  matrix[1][0] = (v[0] * v[1] * (1 - cos_term)) + (v[2] * sin_term);
	  matrix[1][1] = vsq[1] + (vsq[2] + vsq[0])*(cos_term);
	  matrix[1][2] = (v[1] * v[2] * (1 - cos_term)) - (v[0] * sin_term);
	  matrix[2][0] = (v[2] * v[0] * (1 - cos_term)) - (v[1] * sin_term);
	  matrix[2][1] = (v[1] * v[2] * (1 - cos_term)) + (v[0] * sin_term);
	  matrix[2][2] = vsq[2] + (vsq[0] + vsq[1])*(cos_term);
    
	  CVector<> tempvec(0);

	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    for (size_t jDim(0); jDim < NDIM; ++jDim)      
	      tempvec[iDim] += matrix[iDim][jDim] * osret.orientation[jDim];
    
	  osret.orientation = tempvec;
	}
    }
}

Iflt
CLNOrientation::quadraticRootHunter(orientationStreamType LineA, orientationStreamType LineB,
				    Iflt length, Iflt& t_low, Iflt& t_high) const
{
  Iflt working_time = t_low;
  Iflt timescale = 1e-10 
    * length / (length * (LineA.angularVelocity - LineB.angularVelocity).length() 
		+ (LineA.velocity - LineB.velocity).length());
  bool fwdWorking = false;
  
  size_t w = 0;

  while(t_low < t_high)
    {
      //Always try again from the other side
      fwdWorking = !fwdWorking;
    
      if(++w > 1000)
	{
	  I_cerr() << "Window shrunk thousands of times";
	  
	  return working_time;
	}
    
      orientationStreamType A(LineA), B(LineB);
    
      working_time = (fwdWorking ? t_low : t_high);
      performRotation(A, working_time);
      performRotation(B, working_time);
    

      Iflt deltaT;
      {
	Iflt f0 = F_zeroDeriv(A, B),
	  f1 = F_firstDeriv(A, B),
	  halff2 = 0.5 * F_secondDeriv(A, B),
	  halff2max = 0.5 * F_secondDeriv_max(A, B, length);
	
	if (f0 > 0) halff2max = -halff2max;
	
	{
	  Iflt boundEnhancer;
	  // Enhance bound, no point continuing if the bounds are out of bounds
	  if (!quadraticSolution(boundEnhancer, (fwdWorking? ROOT_SMALLEST_POSITIVE : ROOT_SMALLEST_NEGATIVE), 
				 f0, f1, halff2max))
	    break;
	  
	  (fwdWorking ? t_low : t_high) += boundEnhancer;
	}
	
	if (!quadraticSolution(deltaT, ROOT_SMALLEST_POSITIVE, f0, f1, halff2))
	  continue;
      }
      
      if (((working_time + deltaT) > t_high) || ((working_time + deltaT) < t_low))
	continue;
      
      for(size_t i(1000); i != 0; --i)
	{
	  working_time += deltaT;
	  
	  if((working_time > t_high) || (working_time < t_low))
	    break;
	  
	  performRotation(A, deltaT);
	  performRotation(B, deltaT);
	  
	  if (!quadraticSolution(deltaT, ROOT_SMALLEST_EITHER, 
				 F_zeroDeriv(A, B), F_firstDeriv(A, B), 
				 0.5 * F_secondDeriv(A, B)))
	    break;
	  
	  if(fabs(deltaT) <  timescale)
	    return working_time + deltaT;
	}
    }

  return HUGE_VAL;
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

  Iflt factor = std::sqrt(6.0/(length*length));

  CVector<> angVelCrossing;

  for (size_t i = 0; i < Sim->vParticleList.size(); ++i)
    {      
      //Assign the new velocities
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        orientationData[i].orientation[iDim] = Sim->normal_sampler();

      orientationData[i].orientation /= orientationData[i].orientation.length();

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        angVelCrossing[iDim] = Sim->normal_sampler();
      
      orientationData[i].angularVelocity 
        = orientationData[i].orientation.Cross(angVelCrossing).unitVector() 
	* Sim->normal_sampler() * factor;
    }  
}

void 
CLNOrientation::loadParticleXMLData(const XMLNode& XML, std::istream& os)
{
  I_cout() << "Loading Particle Data ";
  fflush(stdout);
  if (XML.isAttributeSet("AttachedBinary")
      && (std::toupper(XML.getAttribute("AttachedBinary")[0]) == 'Y'))
    {
      if (!XML.isAttributeSet("OrientationDataInc")
	  || (std::toupper(XML.getAttribute("OrientationDataInc")[0]) == 'N'))
	D_throw() << "Orientation data is not present in the binary data,"
		  << " cannot load using this liouvillean.";

      Sim->binaryXML = true;
      unsigned long nPart = boost::lexical_cast<unsigned long>(XML.getAttribute("N"));
      boost::progress_display prog(nPart);
      boost::iostreams::filtering_istream base64Convertor;
      
      base64Convertor.push(boost::iostreams::base64_decoder());
      base64Convertor.push(boost::iostreams::base64cleaner_input_filter());
      base64Convertor.push(boost::iostreams::stream_source<std::istream>(os));

      orientationData.resize(nPart);

      for (unsigned long i = 0; i < nPart; ++i)
	{
	  unsigned long ID;
	  CVector<> vel;
	  CVector<> pos;
	  
	  binaryread(base64Convertor, ID);

	  if (i != ID) 
	    D_throw() << "Binary data corruption detected, id's don't match";

	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, vel[iDim]);
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, pos[iDim]);

	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, orientationData[i].orientation[iDim]);

	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, orientationData[i].angularVelocity[iDim]);
	  
	  vel *= Sim->Dynamics.units().unitVelocity();
	  pos *= Sim->Dynamics.units().unitLength();
	  
	  Sim->vParticleList.push_back(CParticle(pos, vel, ID));

	  ++prog;
	}      
    }
  else
    {
      int xml_iter = 0;
      
      unsigned long nPart = XML.nChildNode("Pt");
      boost::progress_display prog(nPart);
      bool outofsequence = false;  
      
      orientationData.resize(nPart);
      
      for (unsigned long i = 0; i < nPart; ++i)
	{
	  XMLNode xBrowseNode = XML.getChildNode("Pt", &xml_iter);
	  
	  if (boost::lexical_cast<unsigned long>
	      (xBrowseNode.getAttribute("ID")) != i)
	    outofsequence = true;
	  
	  CParticle part(xBrowseNode, i);
	  part.scaleVelocity(Sim->Dynamics.units().unitVelocity());
	  part.scalePosition(Sim->Dynamics.units().unitLength());
	  Sim->vParticleList.push_back(part);

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

      if (outofsequence)
	I_cout() << IC_red 
		 << "Particle ID's out of sequence!\n"
		 << IC_red 
		 << "This can result in incorrect capture map loads etc.\n"
		 << IC_red 
		 << "Erase any capture maps in the configuration file so they are regenerated."
		 << IC_reset;            
    }  
}

//! \brief Helper function for writing out data
template<class T>
void 
CLNOrientation::binarywrite(std::ostream& os, const T& val ) const
{
  os.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template<class T>
void 
CLNOrientation::binaryread(std::istream& os, T& val) const
{
  os.read(reinterpret_cast<char*>(&val), sizeof(T));
}

void 
CLNOrientation::outputParticleBin64Data(std::ostream& os) const
{
  if (!Sim->binaryXML)
    return;
  
  
  boost::iostreams::filtering_ostream base64Convertor;
  base64Convertor.push(boost::iostreams::base64_encoder());
  base64Convertor.push(boost::iostreams::line_wrapping_output_filter(80));
  base64Convertor.push(boost::iostreams::stream_sink<std::ostream>(os));
  
  boost::progress_display prog(Sim->lN);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      CParticle tmp(part);
      Sim->Dynamics.BCs().setPBC(tmp.getPosition(), tmp.getVelocity());
      
      tmp.scaleVelocity(1.0 / Sim->Dynamics.units().unitVelocity());
      tmp.scalePosition(1.0 / Sim->Dynamics.units().unitLength());	  

      binarywrite(base64Convertor, tmp.getID());

      for (size_t iDim(0); iDim < NDIM; ++iDim)
	binarywrite(base64Convertor, tmp.getVelocity()[iDim]);

      for (size_t iDim(0); iDim < NDIM; ++iDim)
	binarywrite(base64Convertor, tmp.getPosition()[iDim]);

      for (size_t iDim(0); iDim < NDIM; ++iDim)
	binarywrite(base64Convertor, orientationData[part.getID()].orientation[iDim]);
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	binarywrite(base64Convertor, orientationData[part.getID()].angularVelocity[iDim]);

      ++prog;
    }
}

void 
CLNOrientation::outputParticleXMLData(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("ParticleData")
      << xmlw::attr("N") << Sim->lN
      << xmlw::attr("AttachedBinary") << (Sim->binaryXML ? "Y" : "N")
      << xmlw::attr("OrientationDataInc") << "Y";

  if (!Sim->binaryXML)
    {      
      I_cout() << "Writing Particles ";
      
      boost::progress_display prog(Sim->lN);

      for (unsigned long i = 0; i < Sim->lN; ++i)
	{
	  CParticle tmp(Sim->vParticleList[i]);
	  Sim->Dynamics.BCs().setPBC(tmp.getPosition(), tmp.getVelocity());
	  
	  tmp.scaleVelocity(1.0 / Sim->Dynamics.units().unitVelocity());
	  tmp.scalePosition(1.0 / Sim->Dynamics.units().unitLength());
	  
	  XML << xmlw::tag("Pt") << tmp; 

	  XML << xmlw::tag("O") 
	      << orientationData[i].angularVelocity
	      << xmlw::endtag("O")
	      << xmlw::tag("U")
	      << orientationData[i].orientation
	      << xmlw::endtag("U")
	      << xmlw::endtag("Pt");
	  
	  ++prog;
	}
    }

  XML << xmlw::endtag("ParticleData");
}

size_t 
CLNOrientation::getParticleDOF() const { return NDIM+2; }

Iflt
CLNOrientation::getParticleKineticEnergy(const CParticle& part) const
{
  /***
   * BODGE: Currently taking length as Sim->Dynamics.units().unitLength()
   *        This is *NOT* always  the length
   */
  
  return 0.5 * Sim->Dynamics.getSpecies(part).getMass()
    *(part.getVelocity().square()
      + Sim->Dynamics.units().unitLength() 
      * Sim->Dynamics.units().unitLength() 
      * orientationData[part.getID()].angularVelocity.square() / 12.0);
}
 
void 
CLNOrientation::rescaleSystemKineticEnergy(const Iflt& scale)
{
  Iflt scalefactor(sqrt(scale));

  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    part.getVelocity() *= scalefactor;

  BOOST_FOREACH(rotData& rdat, orientationData)
    rdat.angularVelocity *= scalefactor;
}

CLNOrientation::IfltPair
CLNOrientation::discIntersectionWindow(orientationStreamType A, orientationStreamType B, Iflt length) const
{
  IfltPair retVal;
  
  CVector<> rij = A.position - B.position;
  CVector<> vij = A.velocity - B.velocity;
  
  Sim->Dynamics.BCs().setPBC(rij, vij);
  
  CVector<> Ahat = A.angularVelocity.unitVector();
  Iflt dotproduct = A.angularVelocity.unitVector() % B.angularVelocity.unitVector();
  
  Iflt signChangeTerm = (length/2.0) * sqrt(1.0 - pow(dotproduct, 2.0));
  
  retVal.alpha = ((-1.0 * (rij % Ahat)) - signChangeTerm)/(vij % Ahat);
  retVal.beta  = ((-1.0 * (rij % Ahat)) + signChangeTerm)/(vij % Ahat);
  
  if(retVal.beta < retVal.alpha)
    std::swap(retVal.alpha, retVal.beta);

  return retVal;
}

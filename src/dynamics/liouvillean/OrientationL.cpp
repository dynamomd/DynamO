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
#include "shapes/lines.hpp"

void 
CLNOrientation::initialise() 
{
  CLiouvillean::initialise();

  Iflt sumEnergy(0.0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)  
    sumEnergy += Sim->Dynamics.getSpecies(part).getScalarMomentOfInertia()
    * orientationData[part.getID()].angularVelocity.nrm2();
  
  //Check if any of the species are overridden
  bool hasInertia(false);
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->Dynamics.getSpecies())
    if (dynamic_cast<const CSpecInertia*>(spec.get_ptr()) != NULL)
      hasInertia = true;

  if (!hasInertia)
    D_throw() << "No species have inertia, using the orientational liouvillean is pointless";

  sumEnergy *= 0.5 / Sim->Dynamics.units().unitEnergy();
  
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
  
  CLinesFunc fL(PD.rij, PD.vij,
		orientationData[p1.getID()].angularVelocity,
		orientationData[p2.getID()].angularVelocity,
		orientationData[p1.getID()].orientation,
		orientationData[p2.getID()].orientation);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    //Shift the lower bound up so we don't find the same root again
    t_low += fabs(2.0 * fL.F_firstDeriv())
      / fL.F_secondDeriv_max(length);
  
  //Find window delimited by discs
  std::pair<Iflt,Iflt> dtw = fL.discIntersectionWindow(length);
  
  if(dtw.first > t_low)
    t_low = dtw.first;
  
  if(dtw.second < t_high)
    t_high = dtw.second;
  
  Iflt root = frenkelRootSearch(fL, length, t_low, t_high);

  if (root != HUGE_VAL) 
    { 
      PD.dt = root; 
      return true; 
    }
  else 
    return false;
}


Iflt
CLNOrientation::frenkelRootSearch(const CLinesFunc& fL, Iflt length, 
				  Iflt t_low, Iflt t_high) const
{
  Iflt root = 0.0;
	
  while(t_high > t_low)
    {
      root = quadraticRootHunter(fL, length, t_low, t_high);

      if (root == HUGE_VAL) return HUGE_VAL;
      
      Iflt temp_high = t_high;

      do {
	// Artificial boundary just below root
	CLinesFunc tempfL(fL);
	tempfL.stream(root);

	Iflt Fdoubleprimemax = tempfL.F_secondDeriv_max(length);
	
	temp_high = root - (fabs(2.0 * tempfL.F_firstDeriv())
			    / Fdoubleprimemax);
	
	if ((temp_high < t_low) || (Fdoubleprimemax == 0)) break;
	
	Iflt temp_root = quadraticRootHunter(fL, length, t_low, temp_high);
	
	if (temp_root == HUGE_VAL) 
	  break;
	else 
	    root = temp_root;
	
      } while(temp_high > t_low);
      
      // At this point $root contains earliest valid root guess.
      // Check root validity.
      CLinesFunc tempfL(fL);
      tempfL.stream(root);
      
      std::pair<Iflt,Iflt> cp = tempfL.getCollisionPoints();
      
      if(fabs(cp.first) < length / 2.0 && fabs(cp.second) < length / 2.0)
        return root;
      else
        t_low = root + ((2.0 * fabs(tempfL.F_firstDeriv()))
			/ tempfL.F_secondDeriv_max(length));
    }
    
  return HUGE_VAL;
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

  retVal.rvdot = (retVal.rij | retVal.vijold);

  Iflt KE1before = getParticleKineticEnergy(particle1);
  Iflt KE2before = getParticleKineticEnergy(particle2);

  CLinesFunc fL(retVal.rij, retVal.vijold,
		orientationData[particle1.getID()].angularVelocity,
		orientationData[particle2.getID()].angularVelocity,
		orientationData[particle1.getID()].orientation,
		orientationData[particle2.getID()].orientation);

  Vector uPerp = fL.getu1() ^ fL.getu2();

  uPerp /= uPerp.nrm();

  std::pair<Iflt, Iflt> cp = fL.getCollisionPoints();

  // \Delta {\bf v}_{imp}
  Vector  vr = retVal.vijold
    + (cp.first * fL.getw1() ^ fL.getu1()) 
    - (cp.second * fL.getw2() ^ fL.getu2());
  
  Iflt mass = retVal.particle1_.getSpecies().getMass();
  Iflt inertia = retVal.particle1_.getSpecies().getScalarMomentOfInertia();

  retVal.dP = uPerp
    * (((vr | uPerp) * (1.0 + elasticity))
       / ((2.0 / mass) + ((cp.first * cp.first + cp.second * cp.second) / inertia)));
  
  const_cast<CParticle&>(particle1).getVelocity() -= retVal.dP / mass;
  const_cast<CParticle&>(particle2).getVelocity() += retVal.dP / mass;

  orientationData[particle1.getID()].angularVelocity 
    -= (cp.first / inertia) * (fL.getu1() ^ retVal.dP);

  orientationData[particle2.getID()].angularVelocity 
    += (cp.second / inertia) * (fL.getu2() ^ retVal.dP);

  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->dSysTime;

  return retVal;
}

void
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  part.getPosition() += part.getVelocity() * dt;

  //The Vector copy is required to make sure that the cached
  //orientation doesn't change during calculation
  orientationData[part.getID()].orientation 
    = Rodrigues(orientationData[part.getID()].angularVelocity * dt)
    * Vector(orientationData[part.getID()].orientation);    
}

Iflt
CLNOrientation::quadraticRootHunter(const CLinesFunc& fL, Iflt length, 
				    Iflt& t_low, Iflt& t_high) const
{
  Iflt working_time = t_low;
  Iflt timescale = 1e-10 * length / fL.F_firstDeriv_max(length);
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
    
      working_time = (fwdWorking ? t_low : t_high);
      CLinesFunc tempfL(fL);
      tempfL.stream(working_time);
    
      Iflt deltaT;
      {
	Iflt f0 = tempfL.F_zeroDeriv(),
	  f1 = tempfL.F_firstDeriv(),
	  halff2 = 0.5 * tempfL.F_secondDeriv(),
	  halff2max = 0.5 * tempfL.F_secondDeriv_max(length);
	
	if (f0 > 0) halff2max = -halff2max;
	
	{
	  Iflt boundEnhancer;
	  // Enhance bound, no point continuing if the bounds are out of bounds
	  if (fwdWorking)
	    { if (!quadSolve<ROOT_SMALLEST_POSITIVE>(f0, f1, halff2max, boundEnhancer)) break; }
	  else
	    if (!quadSolve<ROOT_SMALLEST_NEGATIVE>(f0, f1, halff2max, boundEnhancer)) break;
	  
	  (fwdWorking ? t_low : t_high) += boundEnhancer;
	}
	
	if (!quadSolve<ROOT_SMALLEST_POSITIVE>(f0, f1, halff2, deltaT))
	  continue;
      }
      
      if (((working_time + deltaT) > t_high) 
	  || ((working_time + deltaT) < t_low))
	continue;
      
      for(size_t i(1000); i != 0; --i)
	{
	  working_time += deltaT;
	  
	  if((working_time > t_high) || (working_time < t_low))
	    break;
	  
	  tempfL.stream(deltaT);
	  
	  if (!quadSolve<ROOT_SMALLEST_EITHER>(tempfL.F_zeroDeriv(), 
					       tempfL.F_firstDeriv(), 
					       Iflt(0.5 * tempfL.F_secondDeriv()), deltaT))
	    break;
	  
	  if(fabs(deltaT) <  timescale)
	    return working_time + deltaT;
	}
    }

  return HUGE_VAL;
}

C1ParticleData 
CLNOrientation::runAndersenWallCollision(const CParticle& part, 
					 const Vector & vNorm,
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

  Vector  angVelCrossing;

  for (size_t i = 0; i < Sim->vParticleList.size(); ++i)
    {      
      //Assign the new velocities
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        orientationData[i].orientation[iDim] = Sim->normal_sampler();

      orientationData[i].orientation /= orientationData[i].orientation.nrm();

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        angVelCrossing[iDim] = Sim->normal_sampler();
      
      orientationData[i].angularVelocity 
        = orientationData[i].orientation ^ angVelCrossing;
      
      orientationData[i].angularVelocity *= Sim->normal_sampler() * factor 
	/ orientationData[i].angularVelocity.nrm();
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
	  Vector  vel;
	  Vector  pos;
	  
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

	  Iflt oL = orientationData[i].orientation.nrm();
	  
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
  return 0.5 * ((Sim->Dynamics.getSpecies(part).getMass()
    * part.getVelocity().nrm2())
      + (orientationData[part.getID()].angularVelocity.nrm2()
      * Sim->Dynamics.getSpecies(part).getScalarMomentOfInertia()));
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

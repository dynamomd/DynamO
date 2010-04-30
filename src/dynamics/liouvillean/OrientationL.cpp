/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "shapes/frenkelroot.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "shapes/lines.hpp"
#include "../../extcode/binaryHelper.hpp"

void 
LNOrientation::initialise() 
{
  Liouvillean::initialise();

  Iflt sumEnergy(0.0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)  
    sumEnergy += Sim->dynamics.getSpecies(part).getScalarMomentOfInertia()
    * orientationData[part.getID()].angularVelocity.nrm2();
  
  //Check if any of the species are overridden
  bool hasInertia(false);
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->dynamics.getSpecies())
    if (dynamic_cast<const CSpecInertia*>(spec.get_ptr()) != NULL)
      hasInertia = true;

  if (!hasInertia)
    D_throw() << "No species have inertia, using the orientational liouvillean is pointless";

  sumEnergy *= 0.5 / Sim->dynamics.units().unitEnergy();
  
  I_cout() << "System Rotational Energy " << sumEnergy
	   << "\nRotational kT " << sumEnergy / Sim->lN;
}

void
LNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation";
}

bool 
LNOrientation::getLineLineCollision(CPDData& PD, const Iflt& length, 
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

C2ParticleData 
LNOrientation::runLineLineCollision(const CIntEvent& eevent, const Iflt& elasticity, const Iflt& length) const
{
  const CParticle& particle1 = Sim->vParticleList[eevent.getParticle1ID()];
  const CParticle& particle2 = Sim->vParticleList[eevent.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  C2ParticleData retVal(particle1, particle2,
                        Sim->dynamics.getSpecies(particle1),
                        Sim->dynamics.getSpecies(particle2),
                        CORE);
  
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);

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
LNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  part.getPosition() += part.getVelocity() * dt;

  //The Vector copy is required to make sure that the cached
  //orientation doesn't change during calculation
  orientationData[part.getID()].orientation 
    = Rodrigues(orientationData[part.getID()].angularVelocity * dt)
    * Vector(orientationData[part.getID()].orientation);    
}


C1ParticleData 
LNOrientation::runAndersenWallCollision(const CParticle& part, 
					 const Vector & vNorm,
					 const Iflt& sqrtT
					 ) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";
}
  
C1ParticleData 
LNOrientation::randomGaussianEvent(const CParticle& part, 
				    const Iflt& sqrtT) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";  
}

void 
LNOrientation::initLineOrientations(const Iflt& length)
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
LNOrientation::loadParticleXMLData(const XMLNode& XML)
{
  Liouvillean::loadParticleXMLData(XML);

  if (XML.getChildNode("ParticleData").isAttributeSet("AttachedBinary")
      && (std::toupper(XML.getChildNode("ParticleData").getAttribute("AttachedBinary")[0]) == 'Y'))
    {
      boost::progress_display prog(Sim->lN);
      boost::iostreams::filtering_istream base64Convertor;	  
      base64Convertor.push(boost::iostreams::base64_decoder());
      base64Convertor.push(boost::iostreams::base64cleaner_input_filter());
      
      {
	const char* start = XML.getChildNode("AppendedBinaryOrientation").getText();
	base64Convertor.push(boost::make_iterator_range(std::make_pair(start, start + strlen(start))));
      }
      
      orientationData.resize(Sim->lN);
      
      for (unsigned long i = 0; i < Sim->lN; ++i)
	{
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, orientationData[i].orientation[iDim]);
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, orientationData[i].angularVelocity[iDim]);
	  
	  ++prog;
	}
    }
  else
    {
      XMLNode xSubNode = XML.getChildNode("ParticleData");
      int xml_iter = 0;
      
      boost::progress_display prog(Sim->lN);
      
      orientationData.resize(Sim->lN);
      
      for (unsigned long i = 0; i < Sim->lN; ++i)
	{
	  XMLNode xBrowseNode = xSubNode.getChildNode("Pt", &xml_iter);
	  
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
    }
}

void 
LNOrientation::extraXMLParticleData(xmlw::XmlStream& XML, const size_t ID) const
{
  XML << xmlw::tag("O")
      << orientationData[ID].angularVelocity
      << xmlw::endtag("O")
      << xmlw::tag("U")
      << orientationData[ID].orientation
      << xmlw::endtag("U") ;
}

void 
LNOrientation::extraXMLData(xmlw::XmlStream& XML) const
{
  if (Sim->binaryXML)
    {
      XML << xmlw::tag("AppendedBinaryOrientation")
	  << xmlw::chardata();
      
      {
	boost::iostreams::filtering_ostream base64Convertor;
	base64Convertor.push(boost::iostreams::base64_encoder());
	base64Convertor.push(boost::iostreams::line_wrapping_output_filter(80));
	base64Convertor.push(boost::iostreams::stream_sink<std::ostream>(XML.getUnderlyingStream()));
	
	boost::progress_display prog(Sim->lN);
	
	BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
	  {
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, orientationData[part.getID()].orientation[iDim]);
	    
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, orientationData[part.getID()].angularVelocity[iDim]);
	    
	    ++prog;
	  }
      }

      XML << "\n" << xmlw::endtag("AppendedBinaryOrientation");
    }
}

size_t 
LNOrientation::getParticleDOF() const { return NDIM+2; }

Iflt
LNOrientation::getParticleKineticEnergy(const CParticle& part) const
{
  return 0.5 * ((Sim->dynamics.getSpecies(part).getMass()
    * part.getVelocity().nrm2())
      + (orientationData[part.getID()].angularVelocity.nrm2()
      * Sim->dynamics.getSpecies(part).getScalarMomentOfInertia()));
}
 
void 
LNOrientation::rescaleSystemKineticEnergy(const Iflt& scale)
{
  Iflt scalefactor(sqrt(scale));

  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    part.getVelocity() *= scalefactor;

  BOOST_FOREACH(rotData& rdat, orientationData)
    rdat.angularVelocity *= scalefactor;
}


C2ParticleData 
LNOrientation::RoughSpheresColl(const CIntEvent& event, 
				const Iflt& e, 
				const Iflt& et, 
				const Iflt& d2, 
				const EEventType& eType
				) const
{
  const CParticle& particle1 = Sim->vParticleList[event.getParticle1ID()];
  const CParticle& particle2 = Sim->vParticleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  C2ParticleData retVal(particle1, particle2,
			Sim->dynamics.getSpecies(particle1),
			Sim->dynamics.getSpecies(particle2),
			eType);
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
  
  Iflt p1Mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass/(p1Mass+p2Mass);
  
  retVal.rvdot = (retVal.rij | retVal.vijold);

  //The normal impulse
  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());

  Vector eijn = retVal.rij / retVal.rij.nrm();

  //Now the tangential impulse
  Vector gij = retVal.vijold - std::sqrt(d2) * 0.5 
    * ((orientationData[particle1.getID()].angularVelocity
	+ orientationData[particle2.getID()].angularVelocity) ^ eijn);
  
  Vector gijt = (eijn ^ gij) ^ eijn;

  Iflt Jbar = retVal.particle1_.getSpecies().getScalarMomentOfInertia() 
    / (p1Mass * d2 * 0.25);
  
  retVal.dP += (Jbar * (1-et) / (2*(Jbar + 1))) * gijt;

  Iflt KE1before = getParticleKineticEnergy(particle1);
  Iflt KE2before = getParticleKineticEnergy(particle2);

  //This function must edit particles so it overrides the const!
  const_cast<CParticle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<CParticle&>(particle2).getVelocity() += retVal.dP / p2Mass;

  Vector angularVchange = (1-et) / (std::sqrt(d2) * (Jbar+1)) * (eijn ^ gijt);
 
  orientationData[particle1.getID()].angularVelocity
    += angularVchange;
  orientationData[particle2.getID()].angularVelocity 
    += angularVchange;

  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

  return retVal;
}

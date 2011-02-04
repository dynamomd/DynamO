/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../species/inertia.hpp"
#include <magnet/math/matrix.hpp>


void 
LNOrientation::initialise() 
{
  Liouvillean::initialise();

  double sumEnergy(0.0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)  
    sumEnergy += Sim->dynamics.getSpecies(part).getScalarMomentOfInertia()
    * orientationData[part.getID()].angularVelocity.nrm2();
  
  //Check if any of the species are overridden
  bool hasInertia(false);
  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    if (dynamic_cast<const SpInertia*>(spec.get_ptr()) != NULL)
      hasInertia = true;

  if (!hasInertia)
    M_throw() << "No species have inertia, using the orientational liouvillean is pointless";

  sumEnergy *= 0.5 / Sim->dynamics.units().unitEnergy();
  
  I_cout() << "System Rotational Energy " << sumEnergy
	   << "\nRotational kT " << sumEnergy / Sim->N;
}

void
LNOrientation::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type")
      << "NOrientation";
}

bool 
LNOrientation::getLineLineCollision(CPDData& PD, const double& length, 
				     const Particle& p1, const Particle& p2) const
{  
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(p1))
    M_throw() << "Particle1 " << p1.getID() << " is not up to date";

  if (!isUpToDate(p2))
    M_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

  double t_low = 0.0;
  double t_high = PD.dt;
  
  CLinesFunc fL(PD.rij, PD.vij,
		orientationData[p1.getID()].angularVelocity,
		orientationData[p2.getID()].angularVelocity,
		orientationData[p1.getID()].orientation,
		orientationData[p2.getID()].orientation,
		length);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    //Shift the lower bound up so we don't find the same root again
    t_low += fabs(2.0 * fL.F_firstDeriv())
      / fL.F_secondDeriv_max();
  
  //Find window delimited by discs
  std::pair<double,double> dtw = fL.discIntersectionWindow();
  
  if(dtw.first > t_low)
    t_low = dtw.first;
  
  if(dtw.second < t_high)
    t_high = dtw.second;
  
  std::pair<bool,double> root = frenkelRootSearch(fL, t_low, t_high, length * 1e-10);

  if (root.first) 
    { 
      PD.dt = root.second;
      return true; 
    }
  else 
    return false;
}

PairEventData 
LNOrientation::runLineLineCollision(const IntEvent& eevent, const double& elasticity, const double& length) const
{
  const Particle& particle1 = Sim->particleList[eevent.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[eevent.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
                        Sim->dynamics.getSpecies(particle1),
                        Sim->dynamics.getSpecies(particle2),
                        CORE);
  
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);

  retVal.rvdot = (retVal.rij | retVal.vijold);

  double KE1before = getParticleKineticEnergy(particle1);
  double KE2before = getParticleKineticEnergy(particle2);

  CLinesFunc fL(retVal.rij, retVal.vijold,
		orientationData[particle1.getID()].angularVelocity,
		orientationData[particle2.getID()].angularVelocity,
		orientationData[particle1.getID()].orientation,
		orientationData[particle2.getID()].orientation,
		length);

  Vector uPerp = fL.getu1() ^ fL.getu2();

  uPerp /= uPerp.nrm();

  std::pair<double, double> cp = fL.getCollisionPoints();

  // \Delta {\bf v}_{imp}
  Vector  vr = retVal.vijold
    + (cp.first * fL.getw1() ^ fL.getu1()) 
    - (cp.second * fL.getw2() ^ fL.getu2());
  
  double mass = retVal.particle1_.getSpecies().getMass();
  double inertia = retVal.particle1_.getSpecies().getScalarMomentOfInertia();

  retVal.dP = uPerp
    * (((vr | uPerp) * (1.0 + elasticity))
       / ((2.0 / mass) + ((cp.first * cp.first + cp.second * cp.second) / inertia)));
  
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / mass;

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
LNOrientation::streamParticle(Particle& part, const double& dt) const
{
  part.getPosition() += part.getVelocity() * dt;
  
  //The Vector copy is required to make sure that the cached
  //orientation doesn't change during calculation
  orientationData[part.getID()].orientation 
    = Rodrigues(orientationData[part.getID()].angularVelocity * dt)
    * Vector(orientationData[part.getID()].orientation); 
}


ParticleEventData 
LNOrientation::runAndersenWallCollision(const Particle& part, 
					 const Vector & vNorm,
					 const double& sqrtT
					 ) const
{
  M_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";
}
  
ParticleEventData 
LNOrientation::randomGaussianEvent(const Particle& part, 
				    const double& sqrtT) const
{
  M_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";  
}

void 
LNOrientation::initLineOrientations(const double& length)
{
  orientationData.resize(Sim->particleList.size());
  
  I_cout() << "Initialising the line orientations";

  double factor = std::sqrt(6.0/(length*length));

  Vector  angVelCrossing;

  for (size_t i = 0; i < Sim->particleList.size(); ++i)
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
      boost::progress_display prog(Sim->N);
      boost::iostreams::filtering_istream base64Convertor;	  
      base64Convertor.push(boost::iostreams::base64_decoder());
      base64Convertor.push(boost::iostreams::base64cleaner_input_filter());
      
      {
	const char* start = XML.getChildNode("AppendedBinaryOrientation").getText();
	base64Convertor.push(boost::make_iterator_range(std::make_pair(start, start + strlen(start))));
      }
      
      orientationData.resize(Sim->N);
      
      for (unsigned long i = 0; i < Sim->N; ++i)
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
      
      boost::progress_display prog(Sim->N);
      
      orientationData.resize(Sim->N);
      
      for (unsigned long i = 0; i < Sim->N; ++i)
	{
	  XMLNode xBrowseNode = xSubNode.getChildNode("Pt", &xml_iter);
	  
	  orientationData[i].orientation << xBrowseNode.getChildNode("U");
	  orientationData[i].angularVelocity << xBrowseNode.getChildNode("O");
	  
	  double oL = orientationData[i].orientation.nrm();
	  
	  if (!(oL > 0.0))
	    M_throw() << "Particle ID " << i 
		      << " orientation vector is zero!";
	  
	  //Makes the vector a unit vector
	  orientationData[i].orientation /= oL;
	  
	  ++prog;
	}
    }
}

void 
LNOrientation::extraXMLParticleData(xml::XmlStream& XML, const size_t ID) const
{
  XML << xml::tag("O")
      << orientationData[ID].angularVelocity
      << xml::endtag("O")
      << xml::tag("U")
      << orientationData[ID].orientation
      << xml::endtag("U") ;
}

void 
LNOrientation::extraXMLData(xml::XmlStream& XML) const
{
  if (Sim->binaryXML)
    {
      XML << xml::tag("AppendedBinaryOrientation")
	  << xml::chardata();
      
      {
	boost::iostreams::filtering_ostream base64Convertor;
	base64Convertor.push(boost::iostreams::base64_encoder());
	base64Convertor.push(boost::iostreams::line_wrapping_output_filter(80));
	base64Convertor.push(boost::iostreams::stream_sink<std::ostream>(XML.getUnderlyingStream()));
	
	boost::progress_display prog(Sim->N);
	
	BOOST_FOREACH(const Particle& part, Sim->particleList)
	  {
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, orientationData[part.getID()].orientation[iDim]);
	    
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, orientationData[part.getID()].angularVelocity[iDim]);
	    
	    ++prog;
	  }
      }

      XML << "\n" << xml::endtag("AppendedBinaryOrientation");
    }
}

size_t 
LNOrientation::getParticleDOF() const { return NDIM+2; }

double
LNOrientation::getParticleKineticEnergy(const Particle& part) const
{
  return 0.5 * ((Sim->dynamics.getSpecies(part).getMass()
    * part.getVelocity().nrm2())
      + (orientationData[part.getID()].angularVelocity.nrm2()
      * Sim->dynamics.getSpecies(part).getScalarMomentOfInertia()));
}
 
void 
LNOrientation::rescaleSystemKineticEnergy(const double& scale)
{
  double scalefactor(sqrt(scale));

  BOOST_FOREACH(Particle& part, Sim->particleList)
    part.getVelocity() *= scalefactor;

  BOOST_FOREACH(rotData& rdat, orientationData)
    rdat.angularVelocity *= scalefactor;
}


PairEventData 
LNOrientation::RoughSpheresColl(const IntEvent& event, 
				const double& e, 
				const double& et, 
				const double& d2, 
				const EEventType& eType
				) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
			Sim->dynamics.getSpecies(particle1),
			Sim->dynamics.getSpecies(particle2),
			eType);
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
  
  double p1Mass = retVal.particle1_.getSpecies().getMass(); 
  double p2Mass = retVal.particle2_.getSpecies().getMass();
  double mu = p1Mass * p2Mass/(p1Mass+p2Mass);
  
  retVal.rvdot = (retVal.rij | retVal.vijold);

  //The normal impulse
  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());

  Vector eijn = retVal.rij / retVal.rij.nrm();

  //Now the tangential impulse
  Vector gij = retVal.vijold - std::sqrt(d2) * 0.5 
    * ((orientationData[particle1.getID()].angularVelocity
	+ orientationData[particle2.getID()].angularVelocity) ^ eijn);
  
  Vector gijt = (eijn ^ gij) ^ eijn;

  double Jbar = retVal.particle1_.getSpecies().getScalarMomentOfInertia() 
    / (p1Mass * d2 * 0.25);
  
  retVal.dP += (Jbar * (1-et) / (2*(Jbar + 1))) * gijt;

  double KE1before = getParticleKineticEnergy(particle1);
  double KE2before = getParticleKineticEnergy(particle2);

  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;

  Vector angularVchange = (1-et) / (std::sqrt(d2) * (Jbar+1)) * (eijn ^ gijt);
 
  orientationData[particle1.getID()].angularVelocity
    += angularVchange;
  orientationData[particle2.getID()].angularVelocity 
    += angularVchange;

  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

  return retVal;
}

ParticleEventData 
LNOrientation::runRoughWallCollision(const Particle& part, 
				     const Vector & vNorm,
				     const double& e,
				     const double& et,
				     const double& r
				     ) const
{
  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);

  double KE1before = getParticleKineticEnergy(part);

  double p1Mass = retVal.getSpecies().getMass(); 

  double Jbar = retVal.getSpecies().getScalarMomentOfInertia()
    / (p1Mass * r * r);

  Vector gij = part.getVelocity() - r
    * (orientationData[part.getID()].angularVelocity ^ vNorm);
  
  Vector gijt = (vNorm ^ gij) ^ vNorm;
  
  const_cast<Particle&>(part).getVelocity()
    -= (1+e) * (vNorm | part.getVelocity()) * vNorm
    + (Jbar * (1-et) / (Jbar + 1)) * gijt;

  Vector angularVchange = (1-et) / (r * (Jbar+1)) * (vNorm ^ gijt);
 
  orientationData[part.getID()].angularVelocity
    += angularVchange;

  retVal.setDeltaKE(getParticleKineticEnergy(part) - KE1before);
  return retVal; 
}

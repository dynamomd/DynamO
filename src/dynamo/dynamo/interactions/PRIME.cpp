/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/interactions/PRIME.hpp>
#include <dynamo/BC/BC.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/ranges/1RRange.hpp>
#include <dynamo/ranges/2RSingle.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

//This is a anonymous namespace, the contents of the anonymous
//namespace will not be available outside of this file. This is used
//to place all of the settings of the model at the top of this file.
namespace {
  const size_t NH = 0;
  const size_t CH = 1; //Also known as C_\alpha
  const size_t CO = 2; //Also known as C

  //Values taken from
  //http://onlinelibrary.wiley.com/doi/10.1002/prot.1100/full
  //These are the basic bead diameters
  const double _PRIME_diameters[] = {3.3, //NH
				     3.7, //CH
				     4.0};//CO
  
  //This is the fluctuation of the bond distance allowed
  const double _PRIME_bond_tolerance = 0.02;

  //This is a list of the bond lengths in the backbone. These only
  //apply to beads seperated by 1 backbone bond.
  //
  //We make this a symmetric tensor to simplify the lookups. The zero
  //entries should never be used.
  const double _PRIME_bond_lengths[] = 
    /*        NH,    CH,    CO, */
    {/*NH*/0.000, 1.460, 1.330,
     /*CH*/1.460, 0.000, 1.510,
     /*CO*/1.330, 1.510, 0.000}
    ;

  //This is a list of the pseudobond lengths in the backbone. These
  //only apply to beads seperated by TWO backbone bonds.
  //
  //We make this a symmetric tensor to simplify the lookups. The zero
  //entries should never be used.
  const double _PRIME_pseudobond_lengths[] = 
    /*       NH,   CH,   CO, */
    {/*NH*/0.00, 2.41, 2.45,
     /*CH*/2.41, 0.00, 2.45,
     /*CO*/2.45, 2.45, 0.00}
    ;

  //The next two constants are for the interactions between beads
  //separated by THREE backbone bonds.
  //
  //This is the special pseudobond length between the CH-CH groups. It
  //is the only pseudobond at this distance.
  const double _PRIME_CH_CH_pseudobond_length = 3.80;
  //
  //This is the scaling factor of the bead diameters if they are
  //closer than four bonds on the same chain.
  const double _PRIME_near_diameter_scale_factor = 0.75;
}

namespace dynamo {
  size_t 
  IPRIME::getType(const size_t particleID) const
  {
#ifdef DYNAMO_DEBUG
    if ((particleID < startID) || (particleID > endID))
      M_throw() << "Error, the supplied particle ID is out of range";
#endif
    //The types of the ID's are constants defined above.
    //
    //This implementation will need to be improved once we generalise
    //to multiple backbones.
    return (particleID - startID) % 3;
  }

  size_t 
  IPRIME::getDistance(const size_t pID1, const size_t pID2) const
  {
#ifdef DYNAMO_DEBUG
    if ((pID1 < startID) || (pID2 > endID))
      M_throw() << "Error, the supplied particle ID 1 is out of range";

    if ((pID2 < startID) || (pID2 > endID))
      M_throw() << "Error, the supplied particle ID 2 is out of range";
#endif

    //Here, we must be careful to return only positive distances along
    //the chain. Negative values would cause the unsigned size_t to
    //roll over, plus it makes the algorithm easier to implement if
    //this is always positive.
    //
    //This implementation will need to be improved once we generalise
    //to multiple backbones.
    return std::max(pID1, pID2) - std::min(pID1, pID2);
  }


  IPRIME::IPRIME(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ISingleCapture(tmp, NULL) //A temporary value!
  {
    operator<<(XML);
  }

  void IPRIME::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"PRIME"))
      M_throw() << "Attempting to load PRIME from non PRIME entry";
    

    //We don't call this operator, as we custom load our Range (it must be a linear range)
    //Interaction::operator<<(XML);    
    try {
      startID = XML.getAttribute("Start").as<unsigned long>();
      endID = XML.getAttribute("End").as<unsigned long>() + 1;

      range = shared_ptr<C2Range>(new C2RSingle(new RRange(startID, endID)));

      intName = XML.getAttribute("Name");
    }
    catch (boost::bad_lexical_cast &) { M_throw() << "Failed a lexical cast in IPRIME"; }
  }

  Vector
  IPRIME::getGlyphSize(size_t ID, size_t subID) const 
  {
    //Here we return the hard core diameter
    double diam = _PRIME_diameters[getType(ID)];
    return Vector(diam, diam, diam);
  }

  Vector 
  IPRIME::getGlyphPosition(size_t ID, size_t subID) const 
  {
    Vector retval = Sim->particles[ID].getPosition();
    Sim->BCs->applyBC(retval);
    return retval;
  }


  double 
  IPRIME::getExcludedVolume(size_t ID) const 
  {
    //This calculation only includes the volumes which are always
    //excluded (i.e. the hard core)
    double diam = _PRIME_diameters[getType(ID)];
    return diam * diam * diam * M_PI / 6.0; 
  }

  double 
  IPRIME::maxIntDist() const 
  { 
    //This function cannot be implemented until I figure out the bonds
    //and square well distances
    M_throw() << "More work needed";
  }

  void 
  IPRIME::initialise(size_t nID)
  {
    ID = nID;
    //ISingleCapture::initCaptureMap();
  }

  bool IPRIME::captureTest(const Particle& p1, const Particle& p2) const 
  {
    M_throw() << "Not implemented";

//    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
//
//    double d = (_diameter->getProperty(p1.getID())
//		+ _diameter->getProperty(p2.getID())) * 0.5;
//
//    double l = (_lambda->getProperty(p1.getID())
//		+ _lambda->getProperty(p2.getID())) * 0.5;
//
//#ifdef DYNAMO_DEBUG
//    if (Sim->dynamics->sphereOverlap(p1, p2, d))
//      derr << "Warning! Two particles might be overlapping"
//	   << "Overlap is " << Sim->dynamics->sphereOverlap(p1, p2, d) 
//	/ Sim->units.unitLength()
//	   << "\nd = " << d / Sim->units.unitLength() << std::endl;
//#endif
// 
//    return Sim->dynamics->sphereOverlap(p1, p2, l * d);
  }

  IntEvent
  IPRIME::getEvent(const Particle &p1, 
		   const Particle &p2) const 
  {
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    size_t p1Type = getType(p1.getID());
    size_t p2Type = getType(p2.getID());

    //We need to discover what the interaction diameter is for the
    //particles. At first, we assume its equal to the bead diameters
    double diameter = 0.5 * (_PRIME_diameters[p1Type] + _PRIME_diameters[p2Type]);
    bool bonded = false;

    //This treats the special cases if they are 0,1,2, or three backbone
    //bonds apart
    switch (getDistance(p1.getID(), p2.getID()))
      {
      case 0:
	M_throw() << "Invalid backbone distance of 0";
      case 1:
	{//Every type of this interaction is a bonded interaction
	  diameter = _PRIME_bond_lengths[3 * p1Type + p2Type];
	  bonded = true;
	}
	break;
      case 2:
	{//Every type of this interaction is a pseudobond interaction
	  diameter = _PRIME_pseudobond_lengths[3 * p1Type + p2Type];
	  bonded = true;
	}
	break;
      case 3:
	{
	  //Check if this is the special pseudobond
	  if ((p1Type == CH) && (p2Type == CH))
	    {
	      diameter = _PRIME_CH_CH_pseudobond_length;
	      bonded = true;
	    }
	  else
	    //Its a standard bead-bead interaction, but the diameters
	    //are scaled by a factor
	    diameter *= _PRIME_near_diameter_scale_factor;
	}
	break;
      default:
	//The diameter was correctly specified in the start as a
	//bead-bead interaction.
	break;
      }

    //Now that the diameters are known, and we know if they are bonded
    //or not, we can do event detection. If it is bonded
    //(bonded=true), the inner and outer interactions are at:
    //diameter*(1+-_PRIME_bond_tolerance)

    M_throw() << "Not implemented";

//    double d = (_diameter->getProperty(p1.getID())
//		+ _diameter->getProperty(p2.getID())) * 0.5;
//
//    double l = (_lambda->getProperty(p1.getID())
//		+ _lambda->getProperty(p2.getID())) * 0.5;
//
//    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);
//
//    if (isCaptured(p1, p2))
//      {
//	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
//	if (dt != HUGE_VAL)
//	  {
//#ifdef DYNAMO_OverlapTesting
//	if (Sim->dynamics->sphereOverlap(p1, p2, d))
//	  M_throw() << "Overlapping particles found"
//		    << ", particle1 " << p1.getID()
//		    << ", particle2 " << p2.getID()
//		    << "\nOverlap = " 
//		    << Sim->dynamics.getDynamics()
//	    .sphereOverlap(p1, p2, d)
//	    / Sim->units.unitLength();
//#endif	    
//	    retval = IntEvent(p1, p2, dt, CORE, *this);
//	  }
//
//	dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l * d);
//	if (retval.getdt() > dt)
//	    retval = IntEvent(p1, p2, dt, WELL_OUT, *this);
//      }
//    else
//      {
//	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l * d);
//
//      if (dt != HUGE_VAL)
//	{
//#ifdef DYNAMO_OverlapTesting
//	  if (Sim->dynamics->sphereOverlap(p1, p2, l * d))
//	    {
//	      if (Sim->dynamics->sphereOverlap(p1, p2, d))
//		M_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
//			  << "\nparticle1 " << p1.getID() << ", particle2 " 
//			  << p2.getID() << "\nOverlap = " 
//			  << Sim->dynamics->sphereOverlap(p1, p2, d)
//		  / Sim->units.unitLength();
//	      else
//		M_throw() << "Overlapping wells (but not registerd as captured) particles found" 
//			  << "\nparticle1 " << p1.getID() << ", particle2 " 
//			  << p2.getID() << "\nOverlap = " 
//			  << Sim->dynamics->sphereOverlap(p1, p2, l * d)
//		  / Sim->units.unitLength();
//	      
//	    }
//#endif
//	  retval = IntEvent(p1, p2, dt, WELL_IN, *this);
//	}
//      }
//
//    return retval;
  }

  void
  IPRIME::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
  {
    M_throw() << "Not implemented";

//    ++Sim->eventCount;
//
//    double d = (_diameter->getProperty(p1.getID())
//		+ _diameter->getProperty(p2.getID())) * 0.5;
//    double d2 = d * d;
//
//    double e = (_e->getProperty(p1.getID())
//		+ _e->getProperty(p2.getID())) * 0.5;
//
//    double l = (_lambda->getProperty(p1.getID())
//		+ _lambda->getProperty(p2.getID())) * 0.5;
//    double ld2 = d * l * d * l;
//
//    double wd = (_wellDepth->getProperty(p1.getID())
//		 + _wellDepth->getProperty(p2.getID())) * 0.5;
//    switch (iEvent.getType())
//      {
//      case CORE:
//	{
//	  PairEventData retVal(Sim->dynamics->SmoothSpheresColl(iEvent, e, d2, CORE));
//	  Sim->signalParticleUpdate(retVal);
//	
//	  Sim->ptrScheduler->fullUpdate(p1, p2);
//	
//	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
//	    Ptr->eventUpdate(iEvent, retVal);
//
//	  break;
//	}
//      case WELL_IN:
//	{
//	  PairEventData retVal(Sim->dynamics->SphereWellEvent(iEvent, wd, ld2));
//	
//	  if (retVal.getType() != BOUNCE)
//	    addToCaptureMap(p1, p2);
//	
//	  Sim->ptrScheduler->fullUpdate(p1, p2);
//	  Sim->signalParticleUpdate(retVal);
//	
//	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
//	    Ptr->eventUpdate(iEvent, retVal);
//
//
//	  break;
//	}
//      case WELL_OUT:
//	{
//	  PairEventData retVal(Sim->dynamics->SphereWellEvent(iEvent, -wd, ld2));
//	
//	  if (retVal.getType() != BOUNCE)
//	    removeFromCaptureMap(p1, p2);      
//
//	  Sim->signalParticleUpdate(retVal);
//
//	  Sim->ptrScheduler->fullUpdate(p1, p2);
//	
//	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
//	    Ptr->eventUpdate(iEvent, retVal);
//	  break;
//	}
//      default:
//	M_throw() << "Unknown collision type";
//      } 
  }

  void
  IPRIME::checkOverlaps(const Particle& part1, const Particle& part2) const
  {
    M_throw() << "Not implemented";

//    Vector  rij = part1.getPosition() - part2.getPosition();
//    Sim->BCs->applyBC(rij);
//    double r2 = rij.nrm2();
//
//    double d = (_diameter->getProperty(part1.getID())
//		+ _diameter->getProperty(part2.getID())) * 0.5;
//
//    double l = (_lambda->getProperty(part1.getID())
//		+ _lambda->getProperty(part2.getID())) * 0.5;
//  
//    double d2 = d * d;
//    double ld2 = d2 * l * l;
//
//    if (isCaptured(part1, part2)) {
//	if (r2 < d2)
//	  derr << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
//	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
//	       << r2 / pow(Sim->units.unitLength(),2)
//	       << "\nd^2=" 
//	       << d2 / pow(Sim->units.unitLength(),2) << std::endl;
//
//	if (r2 > ld2)
//	  derr << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
//	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
//	       << r2 / pow(Sim->units.unitLength(),2)
//	       << "\n(lambda * d)^2=" 
//	       << ld2 / pow(Sim->units.unitLength(),2) << std::endl;
//      }
//    else 
//      if (r2 < ld2) {
//	  if (r2 < d2)
//	    derr << "Overlap error\n ID1=" << part1.getID() << ", ID2=" << part2.getID() << "\nR_ij^2=" 
//		 << r2 / pow(Sim->units.unitLength(),2) << "\n(d)^2=" << d2 / pow(Sim->units.unitLength(),2) << std::endl;
//	  else
//	    derr << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
//		 << ", ID2=" << part2.getID() << "\nR_ij^2=" << r2 / pow(Sim->units.unitLength(),2)
//		 << "\n(lambda * d)^2=" << ld2 / pow(Sim->units.unitLength(),2) << std::endl;
//	}
  }
  
  void 
  IPRIME::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "PRIME"
    	<< magnet::xml::attr("Name") << intName
	<< magnet::xml::attr("Start") << startID
	<< magnet::xml::attr("End") << endID
      ;
    
    ISingleCapture::outputCaptureMap(XML);  
  }

  double 
  IPRIME::getInternalEnergy() const
  { 
    M_throw() << "Not implemented";

    ////Once the capture maps are loaded just iterate through that determining energies
    //double Energy = 0.0;
    //typedef std::pair<size_t, size_t> locpair;
    //
    //BOOST_FOREACH(const locpair& IDs, captureMap)
    //  Energy += 0.5 * (_wellDepth->getProperty(IDs.first)
    //		       +_wellDepth->getProperty(IDs.second));
    //
    //return -Energy; 
  }

  double IPRIME::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    M_throw() << "Not implemented";
    //return - 0.5 * (_wellDepth->getProperty(p1.getID()) +_wellDepth->getProperty(p2.getID())) * isCaptured(p1, p2);
  }
}

/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/dynamics/interactions/swsequence.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ISWSequence::ISWSequence(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Interaction(tmp, NULL), //A temporary value!
    _unitEnergy(Sim->_properties.getProperty
		(1.0, Property::Units::Energy()))
  { operator<<(XML); }

  void 
  ISWSequence::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "SquareWellSeq"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Lambda") << _lambda->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;

    XML << magnet::xml::tag("Sequence");
  
    for (size_t i = 0; i < sequence.size(); ++i)
      XML << magnet::xml::tag("Element")
	  << magnet::xml::attr("seqID") << i
	  << magnet::xml::attr("Letter") << sequence[i]
	  << magnet::xml::endtag("Element");

    XML << magnet::xml::endtag("Sequence")
	<< magnet::xml::tag("Alphabet");

    for (size_t i = 0; i < alphabet.size(); ++i)
      for (size_t j = i; j < alphabet[i].size(); ++j)
	XML << magnet::xml::tag("Word")
	    << magnet::xml::attr("Letter1") << i
	    << magnet::xml::attr("Letter2") << j
	    << magnet::xml::attr("Depth") << alphabet[i][j] * _unitEnergy->getMaxValue()
	    << magnet::xml::endtag("Word");

    XML << magnet::xml::endtag("Alphabet");

  
    ISingleCapture::outputCaptureMap(XML);  
  }

  void 
  ISWSequence::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"SquareWellSeq"))
      M_throw() << "Attempting to load SquareWell from non SquareWell entry";
  
    Interaction::operator<<(XML);

    try { 
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
      _lambda = Sim->_properties.getProperty(XML.getAttribute("Lambda"),
					     Property::Units::Dimensionless());

      if (XML.hasAttribute("Elasticity"))
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
      else
	_e = Sim->_properties.getProperty(1.0, Property::Units::Dimensionless());

      intName = XML.getAttribute("Name");
      ISingleCapture::loadCaptureMap(XML);

      //Load the sequence
      sequence.clear();
      std::set<size_t> letters;
    
      for (magnet::xml::Node node = XML.getNode("Sequence").fastGetNode("Element");
	   node.valid(); ++node)
	{
	  if (node.getAttribute("seqID").as<size_t>() != sequence.size())
	    M_throw() << "Sequence of letters not in order, missing element " << sequence.size();

	  size_t letter = node.getAttribute("Letter").as<size_t>();
	  letters.insert(letter);
	  sequence.push_back(letter);
	}

      //Initialise all the well depths to 1.0
      alphabet.resize(letters.size());

      BOOST_FOREACH(std::vector<double>& vec, alphabet)
	vec.resize(letters.size(), 0.0);

      for (magnet::xml::Node node = XML.getNode("Alphabet").fastGetNode("Word");
	   node.valid(); ++node)
	{
	  alphabet
	    .at(node.getAttribute("Letter1").as<size_t>())
	    .at(node.getAttribute("Letter2").as<size_t>())
	    = node.getAttribute("Depth").as<double>();

	  alphabet
	    .at(node.getAttribute("Letter2").as<size_t>())
	    .at(node.getAttribute("Letter1").as<size_t>())
	    = node.getAttribute("Depth").as<double>();
	}
    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CISWSequence";
      }
  }

  double 
  ISWSequence::getDiameter(size_t ID, size_t subID) const
  { return _diameter->getProperty(ID); }

  Vector 
  ISWSequence::getPosition(size_t ID, size_t subID) const
  { 
    Vector retval = Sim->particleList[ID].getPosition();
    Sim->dynamics.BCs().applyBC(retval);
    return retval;
  }

  double 
  ISWSequence::getInternalEnergy() const 
  { 
    //Once the capture maps are loaded just iterate through that determining energies
    double Energy = 0.0;
    typedef std::pair<size_t, size_t> locpair;

    BOOST_FOREACH(const locpair& IDs, captureMap)
      Energy += alphabet
      [sequence[IDs.first % sequence.size()]]
      [sequence[IDs.second % sequence.size()]] 
      * 0.5 * (_unitEnergy->getProperty(IDs.first)
	       +_unitEnergy->getProperty(IDs.second));
  
    return -Energy; 
  }

  double 
  ISWSequence::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    return -alphabet[sequence[p1.getID() % sequence.size()]][sequence[p2.getID() % sequence.size()]]
      * 0.5 * (_unitEnergy->getProperty(p1.getID())
	       +_unitEnergy->getProperty(p2.getID()))
      * isCaptured(p1, p2);
  }


  double 
  ISWSequence::getExcludedVolume(size_t ID) const 
  { 
    double diam = _diameter->getProperty(ID);
    return diam * diam * diam * M_PI / 6.0; 
  }


  double 
  ISWSequence::maxIntDist() const 
  { return _diameter->getMaxValue() * _lambda->getMaxValue(); }

  void 
  ISWSequence::initialise(size_t nID)
  {
    ID = nID;
    ISingleCapture::initCaptureMap(Sim->particleList);
  }

  bool 
  ISWSequence::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->dynamics.getInteraction(p1, p2))) != this) return false;

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;
  
#ifdef DYNAMO_DEBUG
    if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d))
      derr << "Warning! Two particles might be overlapping"
	   << "Overlap is " << Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d) 
	/ Sim->dynamics.units().unitLength()
	   << "\nd = " << d / Sim->dynamics.units().unitLength() << std::endl;
#endif
 
    return Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, l * d);
  }

  IntEvent 
  ISWSequence::getEvent(const Particle &p1, 
			const Particle &p2) const 
  {    
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

#ifdef DYNAMO_CollDebug
    std::cerr << "\n Testing p1 = " << p1.getID() << " p2 = " << p2.getID();
#endif
    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    if (isCaptured(p1, p2))
      {
	double dt = Sim->dynamics.getLiouvillean()
	  .SphereSphereInRoot(p1, p2, d);
	if (dt != HUGE_VAL) 
	  {
#ifdef DYNAMO_OverlapTesting
	    //Check that there is no overlap 
	    if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d))
	      M_throw() << "Overlapping particles found" 
			<< ", particle1 " << p1.getID() 
			<< ", particle2 " 
			<< p2.getID() << "\nOverlap = " 
			<< Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d) 
		/ Sim->dynamics.units().unitLength();
#endif	  
	    retval = IntEvent(p1, p2, dt, CORE, *this);
	  }
      
	dt = Sim->dynamics.getLiouvillean().SphereSphereOutRoot(p1, p2, l * d);
	if (retval.getdt() > dt)
	  retval = IntEvent(p1, p2, dt, WELL_OUT, *this);
      }
    else
      {
	double dt = Sim->dynamics.getLiouvillean()
	  .SphereSphereInRoot(p1, p2, l * d);
	if (dt != HUGE_VAL)
	  {
#ifdef DYNAMO_OverlapTesting
	    if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, l * d))
	      {
		if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d))
		  M_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
			    << "\nparticle1 " << p1.getID() << ", particle2 " 
			    << p2.getID() << "\nOverlap = " 
			    << Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d) 
		    / Sim->dynamics.units().unitLength();
		else
		  M_throw() << "Overlapping wells (but not registerd as captured) particles found" 
			    << "\nparticle1 " << p1.getID() << ", particle2 " 
			    << p2.getID() << "\nOverlap = " 
			    << Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, l * d)
		    / Sim->dynamics.units().unitLength();
	  
	      }
#endif
	    retval = IntEvent(p1, p2, dt, WELL_IN, *this);
	  }
      }
    return retval;
  }

  void
  ISWSequence::runEvent(const Particle& p1,
			const Particle& p2,
			const IntEvent& iEvent) const
  {  
    ++Sim->eventCount;

    double e = (_e->getProperty(p1.getID())
		+ _e->getProperty(p2.getID())) * 0.5;

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double d2 = d * d;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;
  
    double ld2 = d * l * d * l;

    switch (iEvent.getType())
      {
      case CORE:
	{
	  PairEventData retVal(Sim->dynamics.getLiouvillean().SmoothSpheresColl(iEvent, e, d2, CORE));
	  Sim->signalParticleUpdate(retVal);
	
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);

	  break;
	}
      case WELL_IN:
	{
	  PairEventData retVal(Sim->dynamics.getLiouvillean()
			       .SphereWellEvent
			       (iEvent, alphabet
				[sequence[p1.getID() % sequence.size()]]
				[sequence[p2.getID() % sequence.size()]] 
				* _unitEnergy->getMaxValue(), 
				ld2));
	
	  if (retVal.getType() != BOUNCE)
	    addToCaptureMap(p1, p2);      

	  Sim->signalParticleUpdate(retVal);

	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);

	  break;
	}
      case WELL_OUT:
	{
	  PairEventData retVal(Sim->dynamics.getLiouvillean()
			       .SphereWellEvent
			       (iEvent, -alphabet
				[sequence[p1.getID() % sequence.size()]]
				[sequence[p2.getID() % sequence.size()]]
				* _unitEnergy->getMaxValue(), 
				ld2));
	
	  if (retVal.getType() != BOUNCE)
	    removeFromCaptureMap(p1, p2);
	
	  Sim->signalParticleUpdate(retVal);

	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);

	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      }
  }

  void
  ISWSequence::checkOverlaps(const Particle& part1, const Particle& part2) const
  {
    Vector  rij = part1.getPosition() - part2.getPosition();
    Sim->dynamics.BCs().applyBC(rij);
    double r2 = rij.nrm2();

    double d = (_diameter->getProperty(part1.getID())
		+ _diameter->getProperty(part2.getID())) * 0.5;

    double d2 = d * d;

    double l = (_lambda->getProperty(part1.getID())
		+ _lambda->getProperty(part2.getID())) * 0.5;
  
    double ld2 = d * l * d * l;

    if (isCaptured(part1, part2))
      {
	if (r2 < d2)
	  derr << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\nd^2=" 
	       << d2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;

	if (r2 > ld2)
	  derr << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\n(lambda * d)^2=" 
	       << ld2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;
      }
    else 
      {
	if (r2 < d2)
	  derr << "Particles overlapping cores without even being captured."
	       << "\nProbably a bad initial configuration."
	       << "\n ID1=" 
	       << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\nd^2=" 
	       << d2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;
	if (r2 < ld2)
	  derr << "Possible missed captured pair in diagnostics\n ID1=" 
	       << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\n(lambda * d)^2=" 
	       << ld2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;
      }
  }
}

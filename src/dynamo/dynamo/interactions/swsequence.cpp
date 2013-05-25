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

#include <dynamo/interactions/swsequence.hpp>
#include <dynamo/BC/BC.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ISWSequence::ISWSequence(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ICapture(tmp, NULL), //A temporary value!
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

  
    ICapture::outputCaptureMap(XML);  
  }

  void 
  ISWSequence::operator<<(const magnet::xml::Node& XML)
  {
    Interaction::operator<<(XML);

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
    ICapture::loadCaptureMap(XML);
    
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
    
    for (std::vector<double>& vec : alphabet)
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

  Vector
  ISWSequence::getGlyphSize(size_t ID, size_t subID) const
  { 
    double diam = _diameter->getProperty(ID);
    return Vector(diam, diam, diam); 
  }

  Vector 
  ISWSequence::getGlyphPosition(size_t ID, size_t subID) const
  { 
    Vector retval = Sim->particles[ID].getPosition();
    Sim->BCs->applyBC(retval);
    return retval;
  }

  double 
  ISWSequence::getInternalEnergy() const 
  { 
    //Once the capture maps are loaded just iterate through that determining energies
    double Energy = 0.0;
    for (const ICapture::value_type& IDs : *this)
      Energy += alphabet
      [sequence[IDs.first.first % sequence.size()]]
      [sequence[IDs.first.second % sequence.size()]] 
      * 0.5 * (_unitEnergy->getProperty(IDs.first.first)
	       +_unitEnergy->getProperty(IDs.first.second));
  
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
    ICapture::initCaptureMap();
  }

  size_t
  ISWSequence::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;
  
#ifdef DYNAMO_DEBUG
    if (Sim->dynamics->sphereOverlap(p1, p2, d))
      derr << "Warning! Two particles might be overlapping"
	   << "Overlap is " << Sim->dynamics->sphereOverlap(p1, p2, d) 
	/ Sim->units.unitLength()
	   << "\nd = " << d / Sim->units.unitLength() << std::endl;
#endif
 
    return Sim->dynamics->sphereOverlap(p1, p2, l * d) > 0;
  }

  IntEvent 
  ISWSequence::getEvent(const Particle &p1, 
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

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    if (isCaptured(p1, p2))
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
	if (dt != HUGE_VAL) 
	  retval = IntEvent(p1, p2, dt, CORE, *this);
      
	dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l * d);
	if (retval.getdt() > dt)
	  retval = IntEvent(p1, p2, dt, STEP_OUT, *this);
      }
    else
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l * d);
	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, STEP_IN, *this);
      }
    return retval;
  }

  void
  ISWSequence::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
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
	  PairEventData retVal(Sim->dynamics->SmoothSpheresColl(iEvent, e, d2, CORE));
	  (*Sim->_sigParticleUpdate)(retVal);
	
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);

	  break;
	}
      case STEP_IN:
	{
	  PairEventData retVal(Sim->dynamics->SphereWellEvent(iEvent, alphabet[sequence[p1.getID() % sequence.size()]][sequence[p2.getID() % sequence.size()]] * _unitEnergy->getMaxValue(), ld2, 1));
	  if (retVal.getType() != BOUNCE) ICapture::add(p1, p2);      
	  (*Sim->_sigParticleUpdate)(retVal);
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	  for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);

	  break;
	}
      case STEP_OUT:
	{
	  PairEventData retVal(Sim->dynamics->SphereWellEvent(iEvent, -alphabet[sequence[p1.getID() % sequence.size()]][sequence[p2.getID() % sequence.size()]]	* _unitEnergy->getMaxValue(), ld2, 0));
	  if (retVal.getType() != BOUNCE) ICapture::remove(p1, p2);
	  (*Sim->_sigParticleUpdate)(retVal);
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	  for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      }
  }

  bool
  ISWSequence::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    double d = (_diameter->getProperty(p1.getID()) + _diameter->getProperty(p2.getID())) * 0.5;
    double l = (_lambda->getProperty(p1.getID()) + _lambda->getProperty(p2.getID())) * 0.5;
    if (isCaptured(p1, p2))
      {
	if (!Sim->dynamics->sphereOverlap(p1, p2, l * d))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " registered as being inside the well at " << l * d / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    
	    return true;
	  }

	if (Sim->dynamics->sphereOverlap(p1, p2, d))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " are inside the well with an inner hard core at " << d / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    
	    return true;
	  }
      }
    else
      if (Sim->dynamics->sphereOverlap(p1, p2, l * d))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " are registered as being outside the well at a distance of " << l * d / Sim->units.unitLength()
		 << " but they are at a distance of "
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << std::endl;
	  
	  return true;
	}
    return false;
  }
}

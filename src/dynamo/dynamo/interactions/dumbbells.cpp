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

#include <dynamo/interactions/dumbbells.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IDumbbells::IDumbbells(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ISingleCapture(tmp, NULL)
  {
    operator<<(XML);
  }

  void 
  IDumbbells::initialise(size_t nID)
  {
    ID = nID; 
    ISingleCapture::initCaptureMap();
  }

  Vector IDumbbells::getGlyphSize(size_t ID, size_t subID) const
  { 
    double l;
    if (subID == 0)
      l = _diamA->getProperty(ID);
    else
      l = _diamB->getProperty(ID);

    return Vector(l,l,l);
  }

  Vector IDumbbells::getGlyphPosition(size_t ID, size_t subID) const
  {
    Vector retval = Sim->particles[ID].getPosition();
    Sim->BCs->applyBC(retval);

    double fracA = _fractionA->getProperty(ID);
    double displacement = ((subID == 0) ? fracA : -(1-fracA)) * _separation->getProperty(ID);
    
    return retval + displacement * Sim->dynamics->getRotData(ID).orientation;
  }


  void 
  IDumbbells::operator<<(const magnet::xml::Node& XML)
  { 
    if (strcmp(XML.getAttribute("Type"),"Lines"))
      M_throw() << "Attempting to load Lines from non Lines entry";
  
    Interaction::operator<<(XML);
  
    try 
      {
	_diamA = Sim->_properties.getProperty(XML.getAttribute("DiameterA"),
					       Property::Units::Length());
	_diamB = Sim->_properties.getProperty(XML.getAttribute("DiameterB"),
					       Property::Units::Length());
	_separation = Sim->_properties.getProperty(XML.getAttribute("Separation"),
						   Property::Units::Length());
	_fractionA = Sim->_properties.getProperty(XML.getAttribute("Fraction"),
						  Property::Units::Dimensionless());
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
	intName = XML.getAttribute("Name");
	ISingleCapture::loadCaptureMap(XML);   
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CIDumbbells";
      }
  }

  double 
  IDumbbells::maxIntDist() const 
  { 
    double l = 2 * std::max(_fractionA->getMaxValue() * _separation->getMaxValue() + _diamA->getMaxValue(), 
			    (1-_fractionA->getMaxValue()) * _separation->getMaxValue() + _diamA->getMaxValue());

    return l; 
  }

  IntEvent 
  IDumbbells::getEvent(const Particle &p1, const Particle &p2) const
  {
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 
  
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  IDumbbells::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
  {
    M_throw() << "Error, cannot run these events yet";
//    PairEventData retval;
//
//    switch (iEvent.getType())
//      {
//      case CORE:
//	{
//	  ++Sim->eventCount;
//	  //We have a line interaction! Run it
//	  double e = (_e->getProperty(p1.getID())
//		      + _e->getProperty(p2.getID())) * 0.5;
//	  double l = (_length->getProperty(p1.getID())
//		      + _length->getProperty(p2.getID())) * 0.5;
//
//	  retval = Sim->dynamics->runLineLineCollision(iEvent, e, l);
//	  break;
//	}
//      case NBHOOD_IN:
//	{
//	  addToCaptureMap(p1, p2);
//	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
//	  iEvent.setType(VIRTUAL);
//	  break;
//	}
//      case NBHOOD_OUT:
//	{
//	  removeFromCaptureMap(p1, p2);
//	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
//	  iEvent.setType(VIRTUAL);
//	  break;
//	}
//      case VIRTUAL:
//	{
//	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
//	  iEvent.setType(VIRTUAL);
//	  break;
//	}
//      default:
//	M_throw() << "Unknown collision type";
//      }
//    
//    Sim->signalParticleUpdate(retval);
//    
//    Sim->ptrScheduler->fullUpdate(p1, p2);
//    
//    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, 
//		  Sim->outputPlugins)
//      Ptr->eventUpdate(iEvent, retval);
  }
   
  void 
  IDumbbells::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Dumbbells"
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("DiameterA") << _diamA->getName()
	<< magnet::xml::attr("DiameterB") << _diamB->getName()
	<< magnet::xml::attr("Separation") << _separation->getName()
	<< magnet::xml::attr("Fraction") << _fractionA->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;

    ISingleCapture::outputCaptureMap(XML);
  }

  bool 
  IDumbbells::captureTest(const Particle& p1, const Particle& p2) const
  {
    return false;
//    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
//
//    double l = (_length->getProperty(p1.getID())
//		+ _length->getProperty(p2.getID())) * 0.5;
//
//    return Sim->dynamics->sphereOverlap(p1, p2, l);
  }

  bool
  IDumbbells::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    return false;
  }
}

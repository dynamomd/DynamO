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

    return retval + ((subID == 0) ? _LA->getProperty(ID) : -_LB->getProperty(ID)) * Sim->dynamics->getRotData(ID).orientation;
  }


  void 
  IDumbbells::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);
  
    try 
      {
	_diamA = Sim->_properties.getProperty(XML.getAttribute("DiameterA"),
					      Property::Units::Length());
	_diamB = Sim->_properties.getProperty(XML.getAttribute("DiameterB"),
					      Property::Units::Length());
	_LA = Sim->_properties.getProperty(XML.getAttribute("LA"),
					   Property::Units::Length());
	_LB = Sim->_properties.getProperty(XML.getAttribute("LB"),
					   Property::Units::Length());
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
    double l = std::max(2 * _LA->getMaxValue() + _diamA->getMaxValue(), 
			2 * _LB->getMaxValue() + _diamB->getMaxValue());
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
  
    const double lA1 = _LA->getProperty(p1.getID()),
      lB1 = _LB->getProperty(p1.getID()),
      diamA1 = _diamA->getProperty(p1.getID()),
      diamB1 = _diamB->getProperty(p1.getID());

    const double lA2 = _LA->getProperty(p2.getID()),
      lB2 = _LB->getProperty(p2.getID()),
      diamA2 = _diamA->getProperty(p2.getID()),
      diamB2 = _diamB->getProperty(p2.getID());

    const double l1 = std::max(lA1 + 0.5 * diamA1, lB1 + 0.5 * diamB1);
    const double l2 = std::max(lA2 + 0.5 * diamA2, lB2 + 0.5 * diamB2);
    const double max_dist = l1 + l2;

    if (isCaptured(p1, p2))
      {
	//Run this to determine when the spheres no longer intersect
	const double upper_limit = Sim->dynamics->SphereSphereOutRoot(p1, p2, max_dist);
	
	//Test the AA pairing
	std::pair<bool, double> current = Sim->dynamics->getOffcentreSpheresCollision(lA1, diamA1, lA2, diamA2,
										      p1, p2, upper_limit, max_dist);
	//Test the BB pairing
	std::pair<bool, double> AB = Sim->dynamics->getOffcentreSpheresCollision(lA1, diamA1, -lB2, diamB2,
										 p1, p2, std::min(upper_limit, current.second), max_dist);
	if (AB.second < current.second)
	  current = AB;

	std::pair<bool, double> BA = Sim->dynamics->getOffcentreSpheresCollision(-lB1, diamB1, lA2, diamA2,
										 p1, p2, std::min(upper_limit, current.second), max_dist);

	if (BA.second < current.second)
	  current = BA;

	std::pair<bool, double> BB = Sim->dynamics->getOffcentreSpheresCollision(-lB1, diamB1, -lB2, diamB2,
										 p1, p2, std::min(upper_limit, current.second), max_dist);

	if (BB.second < current.second)
	  current = BB;

	//Check if they miss each other
	if (current.second == HUGE_VAL)
	  return IntEvent(p1, p2, upper_limit, NBHOOD_OUT, *this);
 
	//Something happens in the time interval
	if (current.first)
	  //Its a collision!
	  return IntEvent(p1, p2, current.second, CORE, *this);
	else
	  //Its a virtual event, we need to recalculate in a bit
	  return IntEvent(p1, p2, current.second, VIRTUAL, *this);
      }
    else 
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l1 + l2);
	if (dt != HUGE_VAL)
	  return IntEvent(p1, p2, dt, NBHOOD_IN, *this);
      }

    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  IDumbbells::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
  {
    PairEventData retval;
    
    switch (iEvent.getType())
      {
      case CORE:
	{
	  ++Sim->eventCount;
	  Sim->dynamics->updateParticlePair(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], CORE);
	  p1.getVelocity() = -p1.getVelocity();
	  p2.getVelocity() = -p2.getVelocity();
	  Sim->dynamics->getRotData(p1).angularVelocity = -Sim->dynamics->getRotData(p1).angularVelocity;
	  Sim->dynamics->getRotData(p2).angularVelocity = -Sim->dynamics->getRotData(p2).angularVelocity;
	  break;
	}
      case NBHOOD_IN:
	{
	  addToCaptureMap(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      case NBHOOD_OUT:
	{
	  removeFromCaptureMap(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      case VIRTUAL:
	{
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      }
    
    Sim->signalParticleUpdate(retval);
    
    Sim->ptrScheduler->fullUpdate(p1, p2);
    
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, 
		  Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, retval);
  }
   
  void 
  IDumbbells::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Dumbbells"
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("DiameterA") << _diamA->getName()
	<< magnet::xml::attr("DiameterB") << _diamB->getName()
	<< magnet::xml::attr("LA") << _LA->getName()
	<< magnet::xml::attr("LB") << _LB->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;

    ISingleCapture::outputCaptureMap(XML);
  }

  bool 
  IDumbbells::captureTest(const Particle& p1, const Particle& p2) const
  {
    return false;
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
    
    const double l1 = std::max(_LA->getProperty(p1.getID()) + 0.5 * _diamA->getProperty(p1.getID()), 
			       _LB->getProperty(p1.getID()) + 0.5 * _diamB->getProperty(p1.getID()));
    const double l2 = std::max(_LA->getProperty(p2.getID()) + 0.5 * _diamA->getProperty(p2.getID()), 
			       _LB->getProperty(p2.getID()) + 0.5 * _diamB->getProperty(p2.getID()));

    return Sim->dynamics->sphereOverlap(p1, p2, l1 + l2);
  }

  namespace {
    inline double overlap(const Vector dist, double diam)
    {
      return std::sqrt(std::max(diam * diam - (dist | dist), 0.0));
    }
  }

  bool
  IDumbbells::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const double lA1 = _LA->getProperty(p1.getID()),
      lB1 = _LB->getProperty(p1.getID()),
      diamA1 = _diamA->getProperty(p1.getID()),
      diamB1 = _diamB->getProperty(p1.getID());
    const Vector director1 = Sim->dynamics->getRotData(ID).orientation;

    const double lA2 = _LA->getProperty(p2.getID()),
      lB2 = _LB->getProperty(p2.getID()),
      diamA2 = _diamA->getProperty(p2.getID()),
      diamB2 = _diamB->getProperty(p2.getID());
    const Vector director2 = Sim->dynamics->getRotData(ID).orientation;

    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);

    if (overlap(r12 + director1 * lA1 - director2 * lA2, (diamA1 + diamA2) / 2)
	|| overlap(r12 + director1 * lA1 + director2 * lB2, (diamA1 + diamB2) / 2)
	|| overlap(r12 - director1 * lB1 - director2 * lA2, (diamB1 + diamA2) / 2)
	|| overlap(r12 - director1 * lB1 + director2 * lB2, (diamB1 + diamB2) / 2))
      {
	if (textoutput)
	  derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
	       << " dumbells are overlapping."
	       << std::endl;
	return true;
      }
    return false;
  }
}

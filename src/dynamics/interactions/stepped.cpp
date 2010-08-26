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

#include "stepped.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../BC/BC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/global.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../2particleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <iomanip>
#include <boost/math/special_functions/pow.hpp>

IStepped::IStepped(DYNAMO::SimData* tmp, 
		     const std::vector<steppair>& vec, C2Range* nR):
  IMultiCapture(tmp,nR),
  steps(vec)
{}

IStepped::IStepped(const XMLNode& XML, DYNAMO::SimData* tmp):
  IMultiCapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

void 
IStepped::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Stepped"))
    D_throw() << "Attempting to load Stepped from non Stepped entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    intName = XML.getAttribute("Name");

    if (XML.nChildNode("Step"))
      {
	int xml_iter = 0;
	long counter = XML.nChildNode("Step");
	for (long i = 0; i < counter; ++i)
	  {
	    XMLNode browseNode = XML.getChildNode("Step",&xml_iter);
	    steps.push_back(steppair(boost::lexical_cast<Iflt>
				     (browseNode.getAttribute("R"))
				     * Sim->dynamics.units().unitLength(),
				     boost::lexical_cast<Iflt>
				     (browseNode.getAttribute("E"))
				     * Sim->dynamics.units().unitEnergy()
				     ));
	    
	  }
      }
    else
      D_throw() << "No steppings defined for stepped potential " 
		<< intName;
    
    std::sort(steps.rbegin(), steps.rend());

    IMultiCapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CIStepped";
    }
}

Interaction* 
IStepped::Clone() const 
{ return new IStepped(*this); }

Iflt 
IStepped::hardCoreDiam() const 
{ return steps.back().first; }

Iflt 
IStepped::maxIntDist() const 
{ return steps.front().first; }

void 
IStepped::rescaleLengths(Iflt scale) 
{ 
  BOOST_FOREACH(steppair& p, steps)
    p.first += scale * p.first;
}

void 
IStepped::initialise(size_t nID)
{
  ID = nID;
  IMultiCapture::initCaptureMap();
  
  runstepdata = steps;

  //Make runstepdata hold r^2 and E_i - E_{i-1}
  BOOST_FOREACH(steppair& pdat, runstepdata)
    pdat.first *= pdat.first;

  Iflt oldE(0.0);
  BOOST_FOREACH(steppair& pdat, runstepdata)
    {
      Iflt tmp(pdat.second);
      pdat.second -= oldE;
      oldE = tmp;
    }

  I_cout() << "Buckets in captureMap " << captureMap.bucket_count()
	   << "\nMax bucket count " << captureMap.max_bucket_count()
	   << "\nload Factor " << captureMap.load_factor()
	   << "\nMax load Factor " << captureMap.max_load_factor();
}

int 
IStepped::captureTest(const Particle& p1, const Particle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  
  Iflt r = rij.nrm();

  for (size_t i(0); i < steps.size(); ++i)
    if (r > steps[i].first) return i;

  return steps.size() - 1;
}

Iflt 
IStepped::getInternalEnergy() const 
{ 
  //Once the capture maps are loaded just iterate through that determining energies
  Iflt Energy = 0.0;

  typedef std::pair<const std::pair<size_t, size_t>, int> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    Energy += steps[IDs.second - 1].second;
  
  return Energy; 
}

IntEvent
IStepped::getEvent(const Particle &p1, 
		    const Particle &p2) const
{
  
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
  
  const_cmap_it capstat = getCMap_it(p1,p2);

  if (capstat == captureMap.end())
    {
      //Not captured, test for capture
      if (Sim->dynamics.getLiouvillean().SphereSphereInRoot
	  (colldat, runstepdata.front().first))
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap
	      (colldat, runstepdata.front().first))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - steps.front().first)
	      /Sim->dynamics.units().unitLength();
#endif
	  
	  return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
	}
    }
  else
    {
      //Within the potential, look for further capture or release
      //First check if there is an inner step to interact with
      //Then check for that event first
      if ((capstat->second < static_cast<int>(runstepdata.size()))
	  && (Sim->dynamics.getLiouvillean().SphereSphereInRoot
	      (colldat, runstepdata[capstat->second].first)))
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap
	      (colldat, runstepdata[capstat->second].first))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - steps[capstat->second].first)
	      /Sim->dynamics.units().unitLength();
#endif
	  
	  return IntEvent(p1, p2, colldat.dt, WELL_IN , *this);
	}
      else if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot
	       (colldat, runstepdata[capstat->second-1].first))
	return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }

  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
IStepped::runEvent(const Particle& p1, 
		    const Particle& p2,
		    const IntEvent& iEvent) const
{
  ++Sim->eventCount;

  switch (iEvent.getType())
    {
    case WELL_OUT:
      {
	cmap_it capstat = getCMap_it(p1,p2);
	
	PairEventData retVal(Sim->dynamics.getLiouvillean().SphereWellEvent
			      (iEvent, runstepdata[capstat->second-1].second, 
			       runstepdata[capstat->second -1].first));
	
	if (retVal.getType() != BOUNCE)
	  if (!(--capstat->second))
	    //capstat is zero so delete
	    captureMap.erase(capstat);

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    case WELL_IN:
      {
	cmap_it capstat = getCMap_it(p1,p2);
	
	if (capstat == captureMap.end())
	  capstat = captureMap.insert
	    (captureMapType::value_type
	     ((p1.getID() < p2.getID())
	      ? cMapKey(p1.getID(), p2.getID())
	      : cMapKey(p2.getID(), p1.getID()),
	      0)).first;
	
	PairEventData retVal = Sim->dynamics.getLiouvillean().SphereWellEvent
	  (iEvent, -runstepdata[capstat->second].second,
	   runstepdata[capstat->second].first);
	
	if (retVal.getType() != BOUNCE)
	  ++(capstat->second);
	else if (!capstat->second)
	  captureMap.erase(capstat);	    
	
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	    
	break;
      }
    default:
      D_throw() << "Unknown collision type";
    } 
}

void
IStepped::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  const_cmap_it capstat = getCMap_it(part1,part2);

  if (captureTest(part1,part2) != capstat->second)
    I_cerr() << "Particle " << part1.getID() << " and Particle " << part2.getID()
	     << "\nFailing as captureTest gives " << captureTest(part1,part2)
	     << "\nAnd recorded value is " << capstat->second;
}
  
void 
IStepped::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Stepped"
      << xml::attr("Name") << intName
      << range;

  BOOST_FOREACH(const steppair& s, steps)
    XML << xml::tag("Step")
	<< xml::attr("R") 
	<< s.first / Sim->dynamics.units().unitLength()
	<< xml::attr("E")
	<< s.second / Sim->dynamics.units().unitEnergy()
	<< xml::endtag("Step");
  
  IMultiCapture::outputCaptureMap(XML);  
}

void 
IStepped::write_povray_desc(const DYNAMO::RGB& rgb, 
				const size_t& specID, 
				std::ostream& os) const
{
  os << "#declare intrep" << ID << "center = " 
     << "sphere {\n <0,0,0> " 
     << steps.back().first * 0.5
     << "\n texture { pigment { color rgb<" << rgb.R << "," << rgb.G 
     << "," << rgb.B << "> }}\nfinish { phong 0.9 phong_size 60 }\n}\n";

  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->particleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "center\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
}

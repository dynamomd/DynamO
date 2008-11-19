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

#include "ghost.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

CSysGhost::CSysGhost(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
  uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
  meanFreeTime(100000),
  Temp(Sim->Dynamics.units().unitEnergy()),
  sqrtTemp(std::sqrt(Sim->Dynamics.units().unitEnergy())),
  tune(false),
  setPoint(0.05),
  eventCount(0),
  setFrequency(50000),
  range(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = GAUSSIAN;
}

CSysGhost::CSysGhost(DYNAMO::SimData* nSim, Iflt mft, Iflt t, 
		     std::string nName):
  CSystem(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
  meanFreeTime(mft),
  Temp(t),
  tune(true),
  setPoint(0.05),
  eventCount(0),
  setFrequency(50000),
  range(new CRAll(Sim))
{
  sysName = nName;
  type = GAUSSIAN;
}

void 
CSysGhost::stream(Iflt ndt)
{
  dt -= ndt;
  if (tune)
    if ((!(Sim->lNColl % setFrequency)) 
	&& (eventCount != 0) && (Sim->lNColl != 0))
      {
	meanFreeTime *= static_cast<Iflt>(eventCount)
	  / (setFrequency * setPoint);
	dt = getGhostt();
	eventCount = 0;
      }
}

void 
CSysGhost::runEvent() const
{
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(locdt))
    D_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->Dynamics.stream(locdt);

  ++eventCount;

  dt = getGhostt();

  unsigned int step = boost::variate_generator
    <DYNAMO::baseRNG&, boost::uniform_int<unsigned int> >
    (Sim->ranGenerator, 
     boost::uniform_int<unsigned int>(0, range->size() - 1))();

  const CParticle& part(Sim->vParticleList[*(range->begin()+step)]);

  //Run the collision and catch the data
  CNParticleData SDat(Sim->Dynamics.Liouvillean().randomGaussianEvent
		      (part, sqrtTemp));
  
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt);
}

void 
CSysGhost::initialise(size_t nID)
{
  ID = nID;
  meanFreeTime /= Sim->lN;
  dt = getGhostt();
  sqrtTemp = sqrt(Temp);
}

void 
CSysGhost::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Andersen"))
    D_throw() << "Attempting to load Andersen from non Andersen entry"; 
  
  try {
    meanFreeTime = boost::lexical_cast<Iflt>(XML.getAttribute("MFT"))
      * Sim->Dynamics.units().unitTime();
    
    Temp = boost::lexical_cast<Iflt>(XML.getAttribute("Temperature")) 
      * Sim->Dynamics.units().unitEnergy();
    
    sysName = XML.getAttribute("Name");

    if (XML.isAttributeSet("SetFrequency") && XML.isAttributeSet("SetPoint"))
      {
	tune = true;
	
	setFrequency = boost::lexical_cast<unsigned long>
	  (XML.getAttribute("SetFrequency"));

	setPoint = boost::lexical_cast<Iflt>(XML.getAttribute("SetPoint"));
      }

    range.set_ptr(CRange::loadClass(XML,Sim));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSysGhost::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "Andersen"
      << xmlw::attr("Name") << sysName
      << xmlw::attr("MFT") << meanFreeTime
    * Sim->lN
    / Sim->Dynamics.units().unitTime()
      << xmlw::attr("Temperature") << Temp 
    / Sim->Dynamics.units().unitEnergy();
  
  if (tune)
    XML << xmlw::attr("SetPoint") << setPoint
	<< xmlw::attr("SetFrequency") << setFrequency;
  
  XML << range
      << xmlw::endtag("System");
}

Iflt 
CSysGhost::getGhostt() const
{ 
  return  - meanFreeTime * log(1.0-uniformRand());
}

Iflt 
CSysGhost::getReducedTemperature() const
{
  return Temp / Sim->Dynamics.units().unitEnergy();
}

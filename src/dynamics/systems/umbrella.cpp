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

#include "umbrella.hpp"
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

CSUmbrella::CSUmbrella(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
  width(1.0),
  w2(width * width),
  range1(NULL),
  range2(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = UMBRELLA;
}

CSUmbrella::CSUmbrella(DYNAMO::SimData* nSim, Iflt nd, 
		       std::string nName, CRange* r1, CRange* r2):
  CSystem(nSim),
  width(nd),
  w2(nd * nd),
  range1(r1),
  range2(r2)
{
  sysName = nName;
  type = UMBRELLA;
}

void 
CSUmbrella::runEvent() const
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

  ++Sim->lNColl;

  CNParticleData SDat;  


  Sim->signalParticleUpdate(SDat);
}

void
CSUmbrella::initialise(size_t nID)
{
  ID = nID;
  
  CPDData partdata(*Sim, *range1, *range2);

  if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(partdata, width))
    dt = partdata.dt;
  else
    dt = HUGE_VAL;
}

void
CSUmbrella::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Umbrella"))
    D_throw() << "Attempting to load Umbrella from a " 
	      << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    sysName = XML.getAttribute("Name");

    width = boost::lexical_cast<Iflt>(XML.getAttribute("Width"))
      * Sim->Dynamics.units().unitLength();

    range1.set_ptr(CRange::loadClass(XML.getChildNode("Range1"), Sim));

    range2.set_ptr(CRange::loadClass(XML.getChildNode("Range2"), Sim));

  }
  catch (boost::bad_lexical_cast &)
    { D_throw() << "Failed a lexical cast in CSUmbrella"; }
}

void 
CSUmbrella::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "Umbrella"
      << xmlw::attr("width") << width / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Name") << sysName
      << xmlw::tag("Range1")
      << range1
      << xmlw::endtag("Range1")
      << xmlw::tag("Range2")
      << range2
      << xmlw::endtag("Range2")
      << xmlw::endtag("System");
}

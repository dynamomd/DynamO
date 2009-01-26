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
  
  /*const C2ParticleData
    SDat(Sim->Dynamics.Liouvillean().DSMCSpheresRun(p1, p2, e, PDat));
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, 0.0);
  
  Sim->signalParticleUpdate(SDat.particle1_);
  Sim->signalParticleUpdate(SDat.particle2_);*/
  
}

void
CSUmbrella::initialise(size_t nID)
{
  ID = nID;
  
  dt = getTimeTillCollide();
}

Iflt
CSUmbrella::getTimeTillCollide() const
{
  CVector<> COMVel1(0), COMVel2(0), COMPos1(0), COMPos2(0);

  Iflt structmass1(0), structmass2(0);

  BOOST_FOREACH(const size_t& ID, *range1)
    {
      structmass1 += 
	Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

      COMVel1 += Sim->vParticleList[ID].getVelocity()
	* Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

      COMPos1 += Sim->vParticleList[ID].getPosition()
	* Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();
    }

  BOOST_FOREACH(const size_t& ID, *range2)
    {
      structmass2 += 
	Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

      COMVel2 += Sim->vParticleList[ID].getVelocity()
	* Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

      COMPos2 += Sim->vParticleList[ID].getPosition()
	* Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();
    }

  COMVel1 /= structmass1;
  COMVel2 /= structmass2;

  COMPos1 /= structmass1;
  COMPos2 /= structmass2;

  CPDData colldat(*Sim, COMPos1 - COMPos2, COMVel1 - COMVel2);

  if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(colldat, w2))
    return colldat.dt;

  return HUGE_VAL;
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

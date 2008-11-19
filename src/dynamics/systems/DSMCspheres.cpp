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

#include "DSMCspheres.hpp"
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

CSDSMCSpheres::CSDSMCSpheres(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
  uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
  range(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSDSMCSpheres::CSDSMCSpheres(DYNAMO::SimData* nSim, Iflt nd, Iflt ntstp, Iflt nChi, 
			     std::string nName):
  CSystem(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
  tstep(ntstp),
  chi(nChi),
  d2(nd * nd),
  diameter(nd),
  range(new CRAll(Sim))
{
  sysName = nName;
  type = DSMC;
}

void 
CSDSMCSpheres::stream(Iflt ndt)
{
  dt -= ndt;
}

void 
CSDSMCSpheres::runEvent() const
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

  dt = tstep;

  BOOST_FOREACH(const CParticle& p1, Sim->vParticleList)
    BOOST_FOREACH(const CParticle& p2, Sim->vParticleList)
    {
      Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
      CPDData colldat(*Sim, p1, p2);


      
      Sim->ptrScheduler->fullUpdate(p1, p2);
    }
  
  //Run the collision and catch the data
  CNParticleData SDat;
  
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt);
}

void 
CSDSMCSpheres::initialise(size_t nID)
{
  ID = nID;
  dt = tstep;
}

void 
CSDSMCSpheres::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"DSMCSpheres"))
    D_throw() << "Attempting to load DSMCSpheres from a " << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = boost::lexical_cast<Iflt>(XML.getAttribute("tStep"))
      * Sim->Dynamics.units().unitTime();
    
    chi = boost::lexical_cast<Iflt>(XML.getAttribute("Chi"));
    
    sysName = XML.getAttribute("Name");

    diameter = boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"))
      * Sim->Dynamics.units().unitLength();

    d2 = diameter * diameter;

    range.set_ptr(CRange::loadClass(XML,Sim));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSDSMCSpheres::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "DSMCSpheres"
      << xmlw::attr("tStep") << tstep / Sim->Dynamics.units().unitTime()
      << xmlw::attr("Chi") << chi
      << xmlw::attr("Diameter") << diameter / Sim->Dynamics.units().unitLength()
      << range
      << xmlw::endtag("System");
}

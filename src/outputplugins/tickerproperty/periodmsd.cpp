/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "periodmsd.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../0partproperty/msd.hpp"
#include "../../extcode/mathtemplates.hpp"
#include <magnet/math/ctime_pow.hpp>
#include "../../dynamics/ranges/1RAll.hpp"

OPPeriodicMSD::OPPeriodicMSD(const DYNAMO::SimData* tmp, const XMLNode&):
  OPTicker(tmp,"PeriodicMSD"),
  TickerCount(0),
  ptrOPMSD(NULL)
{}

void 
OPPeriodicMSD::initialise()
{
  //Load the diffusion class
  ptrOPMSD = Sim->getOutputPlugin<OPMSD>();

  //Now cache a local list of the topology
  BOOST_FOREACH(const magnet::ClonePtr<Topology>& topo, Sim->dynamics.getTopology())
    {
      localpair2 tmpPair;
      tmpPair.first = topo.get_ptr();
      structResults.push_back(tmpPair);
    }
}

void 
OPPeriodicMSD::ticker()
{
  //Don't do this MSD calc often as its expensive! Taking advantage of
  //a modulus of a power of 2 being fast

  if (!(++TickerCount % magnet::math::ctime_pow<2,4>::result))
    {
      results.push_back
	(std::make_pair(Sim->dSysTime / Sim->dynamics.units().unitTime(), 
			ptrOPMSD->calcMSD(CRAll(Sim))));

      BOOST_FOREACH(localpair2& dat, structResults)
	dat.second.push_back(std::make_pair
			     (Sim->dSysTime / Sim->dynamics.units().unitTime(),
			      ptrOPMSD->calcStructMSD(*dat.first)));
    }
}

void
OPPeriodicMSD::output(xml::XmlStream &XML)
{
  XML << xml::tag("PeriodicMSD") 
      << xml::tag("Particle") 
      << xml::chardata();
  
  BOOST_FOREACH(const localpair& myp, results)
    XML << myp.first << " " << myp.second << "\n";
  
  XML << xml::endtag("Particle");
 
  if (!structResults.empty())
    BOOST_FOREACH(localpair2& dat, structResults)
      {
	XML << xml::tag("Structure")	  
	    << xml::attr("Name") <<  dat.first->getName()
	    << xml::chardata();

	BOOST_FOREACH(const localpair& myp, dat.second)
	  XML << myp.first << " " << myp.second << "\n";

	XML << xml::endtag("Structure");	
      }
 
  
  XML << xml::endtag("PeriodicMSD");
}

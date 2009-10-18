/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

COPPeriodicMSD::COPPeriodicMSD(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"PeriodicMSD"),
  TickerCount(0),
  ptrCOPMSD(NULL)
{}

void 
COPPeriodicMSD::initialise()
{
  //Load the diffusion class
  ptrCOPMSD = Sim->getOutputPlugin<COPMSD>();

  //Now cache a local list of the topology
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& topo, Sim->dynamics.getTopology())
    {
      localpair2 tmpPair;
      tmpPair.first = topo.get_ptr();
      structResults.push_back(tmpPair);
    }
}

void 
COPPeriodicMSD::ticker()
{
  //Don't do this MSD calc often as its expensive! Taking advantage of
  //a modulus of a power of 2 being fast

  if (!(++TickerCount % ctime_pow<2,4>::result))
    {
      results.push_back
	(std::make_pair(Sim->dSysTime / Sim->dynamics.units().unitTime(), 
			ptrCOPMSD->calcMSD()));

      BOOST_FOREACH(localpair2& dat, structResults)
	dat.second.push_back(std::make_pair
			     (Sim->dSysTime / Sim->dynamics.units().unitTime(),
			      ptrCOPMSD->calcStructMSD(*dat.first)));
    }
}

void
COPPeriodicMSD::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("PeriodicMSD") 
      << xmlw::tag("Particle") 
      << xmlw::chardata();
  
  BOOST_FOREACH(const localpair& myp, results)
    XML << myp.first << " " << myp.second << "\n";
  
  XML << xmlw::endtag("Particle");
 
  if (!structResults.empty())
    BOOST_FOREACH(localpair2& dat, structResults)
      {
	XML << xmlw::tag("Structure")	  
	    << xmlw::attr("Name") <<  dat.first->getName()
	    << xmlw::chardata();

	BOOST_FOREACH(const localpair& myp, dat.second)
	  XML << myp.first << " " << myp.second << "\n";

	XML << xmlw::endtag("Structure");	
      }
 
  
  XML << xmlw::endtag("PeriodicMSD");
}

/*  dynamo:- Event driven molecular dynamics simulator 
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
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../0partproperty/msd.hpp"
#include "../../dynamics/ranges/1RAll.hpp"
#include <boost/foreach.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/xmlwriter.hpp>

OPPeriodicMSD::OPPeriodicMSD(const dynamo::SimData* tmp, const magnet::xml::Node&):
  OPTicker(tmp,"PeriodicMSD"),
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

  //Species
  speciesData.resize(Sim->dynamics.getSpecies().size());

}

void 
OPPeriodicMSD::ticker()
{
  BOOST_FOREACH(localpair2& dat, structResults)
    dat.second.push_back(std::make_pair
			 (Sim->dSysTime / Sim->dynamics.units().unitTime(),
			  ptrOPMSD->calcStructMSD(*dat.first)));

  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    speciesData[sp->getID()].push_back(std::make_pair(Sim->dSysTime, ptrOPMSD->calcMSD(*(sp->getRange()))));
}

void
OPPeriodicMSD::output(xml::XmlStream &XML)
{
  XML << xml::tag("PeriodicMSD");
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    {
      XML << xml::tag("Species") 
	  << xml::attr("Name") << sp->getName()
	  << xml::chardata();
      
      BOOST_FOREACH(const localpair& dat, speciesData[sp->getID()])
	XML << dat.first / Sim->dynamics.units().unitTime() << " " 
	    << dat.second * 6 << "\n";

      XML << xml::tag("Species");
    }
  
  if (!structResults.empty())
    {
      BOOST_FOREACH(localpair2& dat, structResults)
	{
	  XML << xml::tag("Structure")	  
	      << xml::attr("Name") <<  dat.first->getName()
	      << xml::chardata();
	  
	  BOOST_FOREACH(const localpair& myp, dat.second)
	    XML << myp.first << " " << myp.second << "\n";
	  
	  XML << xml::endtag("Structure");	
	}
    }
 
  
  XML << xml::endtag("PeriodicMSD");
}

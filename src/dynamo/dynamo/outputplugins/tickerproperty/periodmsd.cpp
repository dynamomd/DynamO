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

#include <dynamo/outputplugins/tickerproperty/periodmsd.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simdata.hpp>
#include <dynamo/liouvillean/liouvillean.hpp>
#include <dynamo/outputplugins/0partproperty/msd.hpp>
#include <dynamo/ranges/1RAll.hpp>
#include <boost/foreach.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPPeriodicMSD::OPPeriodicMSD(const dynamo::SimData* tmp, const magnet::xml::Node&):
    OPTicker(tmp,"PeriodicMSD")
  {}

  void 
  OPPeriodicMSD::initialise()
  {
    //Load the diffusion class
    ptrOPMSD = Sim->getOutputPlugin<OPMSD>();

    if (!ptrOPMSD)
      M_throw() << "Periodic MSD plugin requires MSD plugin!";

    //Now cache a local list of the topology
    BOOST_FOREACH(const shared_ptr<Topology>& topo, Sim->topology)
      {
	localpair2 tmpPair;
	tmpPair.first = topo.get();
	structResults.push_back(tmpPair);
      }

    //Species
    speciesData.resize(Sim->species.size());

  }

  void 
  OPPeriodicMSD::ticker()
  {
    BOOST_FOREACH(localpair2& dat, structResults)
      dat.second.push_back(std::make_pair(Sim->dSysTime, ptrOPMSD->calcStructMSD(*dat.first)));

    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      speciesData[sp->getID()].push_back(std::make_pair(Sim->dSysTime, ptrOPMSD->calcMSD(*(sp->getRange()))));
  }

  void
  OPPeriodicMSD::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("PeriodicMSD");
  
    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      {
	XML << magnet::xml::tag("Species") 
	    << magnet::xml::attr("Name") << sp->getName()
	    << magnet::xml::chardata();
      
	BOOST_FOREACH(const localpair& dat, speciesData[sp->getID()])
	  XML << dat.first / Sim->units.unitTime() << " " 
	      << dat.second << "\n";

	XML << magnet::xml::tag("Species");
      }
  
    if (!structResults.empty())
      {
	BOOST_FOREACH(localpair2& dat, structResults)
	  {
	    XML << magnet::xml::tag("Structure")	  
		<< magnet::xml::attr("Name") <<  dat.first->getName()
		<< magnet::xml::chardata();
	  
	    BOOST_FOREACH(const localpair& myp, dat.second)
	      XML << myp.first / Sim->units.unitTime() << " " << myp.second << "\n";
	  
	    XML << magnet::xml::endtag("Structure");	
	  }
      }
 
  
    XML << magnet::xml::endtag("PeriodicMSD");
  }
}

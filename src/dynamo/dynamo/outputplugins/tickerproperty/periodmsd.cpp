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
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/msd.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPPeriodicMSD::OPPeriodicMSD(const dynamo::Simulation* tmp, const magnet::xml::Node&):
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
    for (const shared_ptr<Topology>& topo : Sim->topology)
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
    for (localpair2& dat : structResults)
      dat.second.push_back(std::make_pair(Sim->systemTime, ptrOPMSD->calcStructMSD(*dat.first)));

    for (const shared_ptr<Species>& sp : Sim->species)
      speciesData[sp->getID()].push_back(std::make_pair(Sim->systemTime, ptrOPMSD->calcMSD(*(sp->getRange()))));
  }

  void
  OPPeriodicMSD::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("PeriodicMSD");
  
    for (const shared_ptr<Species>& sp : Sim->species)
      {
	XML << magnet::xml::tag("Species") 
	    << magnet::xml::attr("Name") << sp->getName()
	    << magnet::xml::chardata();
      
	for (const localpair& dat : speciesData[sp->getID()])
	  XML << dat.first / Sim->units.unitTime() << " " 
	      << dat.second / Sim->units.unitArea() << "\n";

	XML << magnet::xml::tag("Species");
      }
  
    if (!structResults.empty())
      {
	for (localpair2& dat : structResults)
	  {
	    XML << magnet::xml::tag("Structure")	  
		<< magnet::xml::attr("Name") <<  dat.first->getName()
		<< magnet::xml::chardata();
	  
	    for (const localpair& myp : dat.second)
	      XML << myp.first / Sim->units.unitTime() << " " << myp.second / Sim->units.unitArea() << "\n";
	  
	    XML << magnet::xml::endtag("Structure");	
	  }
      }
 
  
    XML << magnet::xml::endtag("PeriodicMSD");
  }
}

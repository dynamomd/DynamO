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

#include <dynamo/outputplugins/tickerproperty/chainBondAngles.hpp>
#include <dynamo/include.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/topology/include.hpp>
#include <dynamo/interactions/captures.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <vector>

namespace dynamo {
  OPChainBondAngles::Cdata::Cdata(size_t ID, size_t CL, double bw):
    chainID(ID)
  {
    BondCorrelations.resize(CL-2, magnet::math::Histogram<>(bw));
    BondCorrelationsAvg.resize(CL-2, 0);
    BondCorrelationsSamples.resize(CL-2, 0);
  }

  OPChainBondAngles::OPChainBondAngles(const dynamo::Simulation* tmp, 
				       const magnet::xml::Node& XML):
    OPTicker(tmp,"ChainBondAngles"),
    binwidth(0.0001)
  {
    operator<<(XML);
  }

  void 
  OPChainBondAngles::operator<<(const magnet::xml::Node& XML)
  {
    if (XML.hasAttribute("binwidth"))
      binwidth = XML.getAttribute("binwidth").as<double>();
  }

  void 
  OPChainBondAngles::initialise()
  {
    for(const shared_ptr<Topology>& plugPtr : Sim->topology)
      if (std::dynamic_pointer_cast<TChain>(plugPtr))
	chains.push_back(Cdata(plugPtr->getID(), plugPtr->getMolecules().front()->size(),
			       binwidth));
  }

  void
  OPChainBondAngles::replicaExchange(OutputPlugin& OPPlug)
  {
    std::swap(Sim, static_cast<OPChainBondAngles&>(OPPlug).Sim);
  }

  void 
  OPChainBondAngles::ticker()
  {
    for (Cdata& dat : chains)
      for (const shared_ptr<IDRange>& range : Sim->topology[dat.chainID]->getMolecules())
      if (range->size() > 2)
	{
	  //Walk the polymer
	  for (size_t j = 0; j < range->size()-2; ++j)
	    {
	      Vector  bond1 = Sim->particles[(*range)[j+1]].getPosition()
		- Sim->particles[(*range)[j]].getPosition();

	      bond1 /= bond1.nrm();

	      for (size_t i = j+2; i < range->size(); ++i)
		{
		  Vector  bond2 = Sim->particles[(*range)[i]].getPosition()
		    -Sim->particles[(*range)[i-1]].getPosition();
		
		  bond2 /= bond2.nrm();
		
		  dat.BondCorrelations[i-j-2].addVal(bond1 | bond2);
		  dat.BondCorrelationsAvg[i-j-2] += (bond1 | bond2);
		  ++(dat.BondCorrelationsSamples[i-j-2]);
		}
	    }
	}
  }

  void 
  OPChainBondAngles::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("BondAngleCorrelators");
  
    for (Cdata& dat : chains)
      {
	XML << magnet::xml::tag("Chain")
	    << magnet::xml::attr("Name") << Sim->topology[dat.chainID]->getName();
            
	size_t Nc = Sim->topology[dat.chainID]
	  ->getMolecules().front()->size() - 2;
      
	for (size_t i = 0; i < Nc; ++i)
	  {
	    XML << magnet::xml::tag("Hist") << magnet::xml::attr("Avg") 
		<< (dat.BondCorrelationsAvg[i] / dat.BondCorrelationsSamples[i]);

	    dat.BondCorrelations[i].outputHistogram(XML, 1.0);

	    XML << magnet::xml::endtag("Hist");
	  }
      
	XML << magnet::xml::endtag("Chain");
      }

    XML << magnet::xml::endtag("BondAngleCorrelators");
  }
}

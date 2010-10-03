/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "chainBondLength.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/ranges/1range.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/topology/include.hpp"
#include "../../dynamics/interactions/captures.hpp"

OPChainBondLength::Cdata::Cdata(size_t ID, size_t CL):
  chainID(ID)
{
  BondLengths.resize(CL-1, C1DHistogram(0.0001));
}

OPChainBondLength::OPChainBondLength(const DYNAMO::SimData* tmp, const XMLNode&):
  OPTicker(tmp,"ChainBondLength")
{}

void 
OPChainBondLength::initialise()
{
  BOOST_FOREACH(const magnet::ClonePtr<Topology>& plugPtr, 
		Sim->dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(Cdata(plugPtr->getID(), 
			     plugPtr->getMolecules().front()->size()));
}

void 
OPChainBondLength::changeSystem(OutputPlugin* OPPlug)
{
  std::swap(Sim, static_cast<OPChainBondLength*>(OPPlug)->Sim);
}

void 
OPChainBondLength::ticker()
{
  BOOST_FOREACH(Cdata& dat, chains)
    BOOST_FOREACH(const magnet::ClonePtr<CRange>& range, 
		  Sim->dynamics.getTopology()[dat.chainID]->getMolecules())
    if (range->size() > 2)
      //Walk the polymer
      for (size_t j = 0; j < range->size()-1; ++j)
	dat.BondLengths[j].addVal
	  ((Sim->particleList[(*range)[j+1]].getPosition()
	    - Sim->particleList[(*range)[j]].getPosition()).nrm());
}

void 
OPChainBondLength::output(xml::XmlStream& XML)
{
  XML << xml::tag("BondAngleLength");
  
  BOOST_FOREACH(Cdata& dat, chains)
    {
      XML << xml::tag("Chain")
	  << xml::attr("Name") 
	  << Sim->dynamics.getTopology()[dat.chainID]->getName();
            
      size_t Nc = Sim->dynamics.getTopology()[dat.chainID]
	->getMolecules().front()->size() - 1;
      
      for (size_t i = 0; i < Nc; ++i)
	dat.BondLengths[i].outputHistogram(XML, 1.0/Sim->dynamics.units().unitLength());
      
      XML << xml::endtag("Chain");
    }
  
  XML << xml::endtag("BondAngleLength");
}

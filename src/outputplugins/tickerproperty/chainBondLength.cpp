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

#include "chainBondLength.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../dynamics/ranges/1range.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/topology/include.hpp"
#include "../../dynamics/interactions/captures.hpp"

COPChainBondLength::Cdata::Cdata(size_t ID, size_t CL):
  chainID(ID)
{
  BondLengths.resize(CL-1, C1DHistogram(0.0001));
}

COPChainBondLength::COPChainBondLength(const DYNAMO::SimData* tmp):
  COPTicker(tmp,"ChainBondLength")
{}

void 
COPChainBondLength::initialise()
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, 
		Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(Cdata(plugPtr->getID(), 
			     plugPtr->getMolecules().front()->size()));
}

void 
COPChainBondLength::changeSystem(COutputPlugin* COPPlug)
{
  std::swap(Sim, static_cast<COPChainBondLength*>(COPPlug)->Sim);
}

void 
COPChainBondLength::ticker()
{
  BOOST_FOREACH(Cdata& dat, chains)
    BOOST_FOREACH(const smrtPlugPtr<CRange>& range, 
		  Sim->Dynamics.getTopology()[dat.chainID]->getMolecules())
    if (range->size() > 2)
      {
	//Update the particles
	BOOST_FOREACH(const unsigned long& id, *range)
	  Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
	
	//Walk the polymer
	for (size_t j = 0; j < range->size()-1; ++j)
	  dat.BondLengths[j].addVal
	    ((Sim->vParticleList[(*range)[j+1]].getPosition()
	      - Sim->vParticleList[(*range)[j]].getPosition()).length());
      }
}

void 
COPChainBondLength::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("BondAngleLength");
  
  BOOST_FOREACH(Cdata& dat, chains)
    {
      XML << xmlw::tag("Chain")
	  << xmlw::attr("Name") 
	  << Sim->Dynamics.getTopology()[dat.chainID]->getName();
            
      size_t Nc = Sim->Dynamics.getTopology()[dat.chainID]
	->getMolecules().front()->size() - 1;
      
      for (size_t i = 0; i < Nc; ++i)
	dat.BondLengths[i].outputHistogram(XML, 1.0/Sim->Dynamics.units().unitLength());
      
      XML << xmlw::endtag("Chain");
    }
  
  XML << xmlw::endtag("BondAngleLength");
}

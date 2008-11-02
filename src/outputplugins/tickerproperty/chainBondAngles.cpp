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

#include "chainBondAngles.hpp"
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

COPChainBondAngles::Cdata::Cdata(size_t ID, size_t CL, Iflt bw):
  chainID(ID)
{
  BondCorrelations.resize(CL-2, C1DHistogram(bw));
  BondCorrelationsAvg.resize(CL-2, 0);
  BondCorrelationsSamples.resize(CL-2, 0);
}

COPChainBondAngles::COPChainBondAngles(const DYNAMO::SimData* tmp, 
				       const XMLNode& XML):
  COPTicker(tmp,"ChainBondAngles"),
  binwidth(0.0001)
{
  operator<<(XML);
}

void 
COPChainBondAngles::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth"))
	binwidth = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPChainBondAngles";
    }
}

void 
COPChainBondAngles::initialise()
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, 
		Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(Cdata(plugPtr->getID(), 
			     plugPtr->getMolecules().front()->size(),
			     binwidth));
}

void 
COPChainBondAngles::changeSystem(COutputPlugin* COPPlug)
{
  std::swap(Sim, static_cast<COPChainBondAngles*>(COPPlug)->Sim);
}

void 
COPChainBondAngles::ticker()
{
  BOOST_FOREACH(Cdata& dat,chains)
    BOOST_FOREACH(const smrtPlugPtr<CRange>& range, 
		  Sim->Dynamics.getTopology()[dat.chainID]->getMolecules())
    if (range->size() > 2)
      {
	//Update the particles
	BOOST_FOREACH(const unsigned long& id, *range)
	  Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
	
	//Walk the polymer
	for (size_t j = 0; j < range->size()-2; ++j)
	  {
	    CVector<> bond1 = Sim->vParticleList[(*range)[j+1]].getPosition()
	      - Sim->vParticleList[(*range)[j]].getPosition();

	    bond1 /= bond1.length();

	    for (size_t i = j+2; i < range->size(); ++i)
	      {
		CVector<> bond2 = Sim->vParticleList[(*range)[i]].getPosition()
		  -Sim->vParticleList[(*range)[i-1]].getPosition();
		
		bond2 /= bond2.length();
		
		dat.BondCorrelations[i-j-2].addVal(bond1 % bond2);
		dat.BondCorrelationsAvg[i-j-2] += (bond1 % bond2);
		++(dat.BondCorrelationsSamples[i-j-2]);
	      }
	  }
      }
}

void 
COPChainBondAngles::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("BondAngleCorrelators");
  
  BOOST_FOREACH(Cdata& dat, chains)
    {
      XML << xmlw::tag("Chain")
	  << xmlw::attr("Name") << Sim->Dynamics.getTopology()[dat.chainID]->getName();
            
      size_t Nc = Sim->Dynamics.getTopology()[dat.chainID]
	->getMolecules().front()->size() - 2;
      
      for (size_t i = 0; i < Nc; ++i)
	{
	  XML << xmlw::tag("Hist") << xmlw::attr("Avg") 
	      << (dat.BondCorrelationsAvg[i] / dat.BondCorrelationsSamples[i]);

	  dat.BondCorrelations[i].outputHistogram(XML, 1.0);

	  XML << xmlw::endtag("Hist");
	}
      
      XML << xmlw::endtag("Chain");
    }

  XML << xmlw::endtag("BondAngleCorrelators");
}

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

OPChainBondAngles::Cdata::Cdata(size_t ID, size_t CL, Iflt bw):
  chainID(ID)
{
  BondCorrelations.resize(CL-2, C1DHistogram(bw));
  BondCorrelationsAvg.resize(CL-2, 0);
  BondCorrelationsSamples.resize(CL-2, 0);
}

OPChainBondAngles::OPChainBondAngles(const DYNAMO::SimData* tmp, 
				       const XMLNode& XML):
  OPTicker(tmp,"ChainBondAngles"),
  binwidth(0.0001)
{
  operator<<(XML);
}

void 
OPChainBondAngles::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth"))
	binwidth = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPChainBondAngles";
    }
}

void 
OPChainBondAngles::initialise()
{
  BOOST_FOREACH(const ClonePtr<Topology>& plugPtr, 
		Sim->dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(Cdata(plugPtr->getID(), 
			     plugPtr->getMolecules().front()->size(),
			     binwidth));
}

void 
OPChainBondAngles::changeSystem(OutputPlugin* OPPlug)
{
  std::swap(Sim, static_cast<OPChainBondAngles*>(OPPlug)->Sim);
}

void 
OPChainBondAngles::ticker()
{
  BOOST_FOREACH(Cdata& dat,chains)
    BOOST_FOREACH(const ClonePtr<CRange>& range, 
		  Sim->dynamics.getTopology()[dat.chainID]->getMolecules())
    if (range->size() > 2)
      {
	//Walk the polymer
	for (size_t j = 0; j < range->size()-2; ++j)
	  {
	    Vector  bond1 = Sim->particleList[(*range)[j+1]].getPosition()
	      - Sim->particleList[(*range)[j]].getPosition();

	    bond1 /= bond1.nrm();

	    for (size_t i = j+2; i < range->size(); ++i)
	      {
		Vector  bond2 = Sim->particleList[(*range)[i]].getPosition()
		  -Sim->particleList[(*range)[i-1]].getPosition();
		
		bond2 /= bond2.nrm();
		
		dat.BondCorrelations[i-j-2].addVal(bond1 | bond2);
		dat.BondCorrelationsAvg[i-j-2] += (bond1 | bond2);
		++(dat.BondCorrelationsSamples[i-j-2]);
	      }
	  }
      }
}

void 
OPChainBondAngles::output(xml::XmlStream& XML)
{
  XML << xml::tag("BondAngleCorrelators");
  
  BOOST_FOREACH(Cdata& dat, chains)
    {
      XML << xml::tag("Chain")
	  << xml::attr("Name") << Sim->dynamics.getTopology()[dat.chainID]->getName();
            
      size_t Nc = Sim->dynamics.getTopology()[dat.chainID]
	->getMolecules().front()->size() - 2;
      
      for (size_t i = 0; i < Nc; ++i)
	{
	  XML << xml::tag("Hist") << xml::attr("Avg") 
	      << (dat.BondCorrelationsAvg[i] / dat.BondCorrelationsSamples[i]);

	  dat.BondCorrelations[i].outputHistogram(XML, 1.0);

	  XML << xml::endtag("Hist");
	}
      
      XML << xml::endtag("Chain");
    }

  XML << xml::endtag("BondAngleCorrelators");
}

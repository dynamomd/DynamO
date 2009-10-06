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

#include "chainContactMap.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../dynamics/ranges/1range.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/topology/include.hpp"
#include "../../dynamics/interactions/captures.hpp"

COPCContactMap::Cdata::Cdata(const CTChain* ptr, unsigned long nMolRange):
  chainPtr(ptr),
  array(new unsigned long[nMolRange * nMolRange]),
  counter(0), chainlength(nMolRange)
{
  for (unsigned long i = 0; i < nMolRange * nMolRange; i++)
    array[i] = 0;
}

COPCContactMap::COPCContactMap(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"ContactMap")
{}

void 
COPCContactMap::initialise()
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(Cdata(dynamic_cast<const CTChain*>(plugPtr.get_ptr()), 
			     plugPtr->getMolecules().front()->size()));
}

void 
COPCContactMap::changeSystem(COutputPlugin* COPPlug)
{
  std::swap(Sim, static_cast<COPCContactMap*>(COPPlug)->Sim);

  BOOST_FOREACH(Cdata& dat, chains)
    {
      try {
	const CTopology* tmpPtr = Sim->Dynamics.getTopology(dat.chainPtr->getName()).get_ptr();
	dat.chainPtr = dynamic_cast<const CTChain*>(tmpPtr);
      } catch (std::exception&)
	{
	  D_throw() << "On changing the system COPCContactMap could not find the topology \"" 
		    << dat.chainPtr->getName() << "\"\n in the new system";
	}
      
      if (dat.chainPtr == NULL)
	D_throw() << "On changing the system COPCContactMap found the topology but failed to upcast!";
    }
}

void 
COPCContactMap::ticker()
{
  BOOST_FOREACH(Cdata& dat,chains)
    BOOST_FOREACH(const smrtPlugPtr<CRange>& range,  dat.chainPtr->getMolecules())
    {
      dat.counter++;
      for (unsigned long i = 0; i < dat.chainlength; i++)
	{
	  const CParticle& part1 = Sim->vParticleList[(*range)[i]];
	 
	  for (unsigned long j = i+1; j < dat.chainlength; j++)
	    {
	      const CParticle& part2 = Sim->vParticleList[(*range)[j]];

	      BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, Sim->Dynamics.getInteractions())
		if (ptr->isInteraction(part1,part2))
		  if (dynamic_cast<const CICapture*>(ptr.get_ptr()) != NULL)
		    if (static_cast<const CICapture*>(ptr.get_ptr())
			->isCaptured(part1,part2))
		      dat.array[i * dat.chainlength + j]++;
	    }
	}
    }
}

void 
COPCContactMap::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("ContactMap");
  
  BOOST_FOREACH(Cdata& dat, chains)
    {
      //Copying
      
      for (unsigned long i = 0; i < dat.chainlength; i++)
	{
	  for (unsigned long j = i+1; j < dat.chainlength; j++)
	    dat.array[j * dat.chainlength + i] = dat.array[i * dat.chainlength + j];
	}

      XML << xmlw::tag(dat.chainPtr->getName().c_str())
	  << xmlw::chardata();

      //Have to draw the boxes so it renders correctly, hence the doubling up      
      for (unsigned long i = 0; i < dat.chainlength; i++)
	{
	  for (unsigned long j = 0; j < dat.chainlength; j++)
	    XML << i - 0.5 << " " << j - 0.5 << " " << (((Iflt) dat.array[i * dat.chainlength + j])/((Iflt) dat.counter)) << "\n"
		<< i - 0.5 << " " << j + 0.5 << " " << (((Iflt) dat.array[i * dat.chainlength + j])/((Iflt) dat.counter)) << "\n";
	  XML << "\n";

	  for (unsigned long j = 0; j < dat.chainlength; j++)
	    XML << i + 0.5 << " " << j - 0.5 << " " << (((Iflt) dat.array[i * dat.chainlength + j])/((Iflt) dat.counter)) << "\n"
		<< i + 0.5 << " " << j + 0.5 << " " << (((Iflt) dat.array[i * dat.chainlength + j])/((Iflt) dat.counter)) << "\n";
	  XML << "\n";

	}
      
            
      XML << xmlw::endtag(dat.chainPtr->getName().c_str());
    }

  XML << xmlw::endtag("ContactMap");
}

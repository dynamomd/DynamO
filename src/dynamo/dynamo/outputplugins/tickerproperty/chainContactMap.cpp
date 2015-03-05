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

#include <dynamo/outputplugins/tickerproperty/chainContactMap.hpp>
#include <dynamo/include.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <dynamo/interactions/captures.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
  OPCContactMap::Cdata::Cdata(const TChain* ptr, unsigned long nMolRange):
    chainPtr(ptr),
    array(new unsigned long[nMolRange * nMolRange]),
    counter(0), chainlength(nMolRange)
  {
    for (unsigned long i = 0; i < nMolRange * nMolRange; i++)
      array[i] = 0;
  }

  OPCContactMap::OPCContactMap(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OPTicker(tmp,"ContactMap")
  {}

  void 
  OPCContactMap::initialise()
  {
    for (const shared_ptr<Topology>& plugPtr : Sim->topology)
      if (std::dynamic_pointer_cast<TChain>(plugPtr))
	chains.push_back(Cdata(static_cast<const TChain*>(plugPtr.get()), plugPtr->getMolecules().front()->size()));
  }

  void 
  OPCContactMap::replicaExchange(OutputPlugin& OPPlug)
  {
    auto& op = static_cast<OPCContactMap&>(OPPlug);
    
#ifdef DYNAMO_DEBUG    
    if (chains.size() != op.chains.size())
      M_throw() << "Chain mismatch in replica exchange";
#endif

    for (size_t i(0); i < chains.size(); ++i)
      {
	std::swap(chains[i].array, op.chains[i].array);
	std::swap(chains[i].counter, op.chains[i].counter);
      }
  }

  void 
  OPCContactMap::ticker()
  {
    for (Cdata& dat : chains)
      for (const shared_ptr<IDRange>& range : dat.chainPtr->getMolecules())
      {
	dat.counter++;
	for (unsigned long i = 0; i < dat.chainlength; i++)
	  {
	    const Particle& part1 = Sim->particles[(*range)[i]];
	 
	    for (unsigned long j = i+1; j < dat.chainlength; j++)
	      {
		const Particle& part2 = Sim->particles[(*range)[j]];

		for (const shared_ptr<Interaction>& ptr : Sim->interactions)
		  if (ptr->isInteraction(part1,part2))
		    if (std::dynamic_pointer_cast<ICapture>(ptr))
		      if (dynamic_cast<const ICapture*>(ptr.get())->isCaptured(part1,part2))
			dat.array[i * dat.chainlength + j]++;
	      }
	  }
      }
  }

  void 
  OPCContactMap::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("ContactMap");
  
    for (Cdata& dat : chains)
      {
	//Copying
      
	for (unsigned long i = 0; i < dat.chainlength; i++)
	  {
	    for (unsigned long j = i+1; j < dat.chainlength; j++)
	      dat.array[j * dat.chainlength + i] = dat.array[i * dat.chainlength + j];
	  }

	XML << magnet::xml::tag(dat.chainPtr->getName().c_str())
	    << magnet::xml::chardata();

	for (size_t i = 0; i < dat.chainlength; i++)
	  {
	    for (size_t j = 0; j < dat.chainlength; ++j)
	      XML << double(dat.array[i * dat.chainlength + j]) / double(dat.counter) << " ";
	    XML << "\n";
	  }
	XML << magnet::xml::endtag(dat.chainPtr->getName().c_str());
      }

    XML << magnet::xml::endtag("ContactMap");
  }
}

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

#include <dynamo/outputplugins/tickerproperty/chainContactMap.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/topology/include.hpp>
#include <dynamo/dynamics/interactions/captures.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
  OPCContactMap::Cdata::Cdata(const CTChain* ptr, unsigned long nMolRange):
    chainPtr(ptr),
    array(new unsigned long[nMolRange * nMolRange]),
    counter(0), chainlength(nMolRange)
  {
    for (unsigned long i = 0; i < nMolRange * nMolRange; i++)
      array[i] = 0;
  }

  OPCContactMap::OPCContactMap(const dynamo::SimData* tmp, const magnet::xml::Node&):
    OPTicker(tmp,"ContactMap")
  {}

  void 
  OPCContactMap::initialise()
  {
    BOOST_FOREACH(const magnet::ClonePtr<Topology>& plugPtr, Sim->dynamics.getTopology())
      if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
	chains.push_back(Cdata(dynamic_cast<const CTChain*>(plugPtr.get_ptr()), 
			       plugPtr->getMolecules().front()->size()));
  }

  void 
  OPCContactMap::changeSystem(OutputPlugin* OPPlug)
  {
    std::swap(Sim, static_cast<OPCContactMap*>(OPPlug)->Sim);

    BOOST_FOREACH(Cdata& dat, chains)
      {
	try {
	  const Topology* tmpPtr = Sim->dynamics.getTopology(dat.chainPtr->getName()).get_ptr();
	  dat.chainPtr = dynamic_cast<const CTChain*>(tmpPtr);
	} catch (std::exception&)
	  {
	    M_throw() << "On changing the system OPCContactMap could not find the topology \"" 
		      << dat.chainPtr->getName() << "\"\n in the new system";
	  }
      
	if (dat.chainPtr == NULL)
	  M_throw() << "On changing the system OPCContactMap found the topology but failed to upcast!";
      }
  }

  void 
  OPCContactMap::ticker()
  {
    BOOST_FOREACH(Cdata& dat,chains)
      BOOST_FOREACH(const magnet::ClonePtr<CRange>& range,  dat.chainPtr->getMolecules())
      {
	dat.counter++;
	for (unsigned long i = 0; i < dat.chainlength; i++)
	  {
	    const Particle& part1 = Sim->particleList[(*range)[i]];
	 
	    for (unsigned long j = i+1; j < dat.chainlength; j++)
	      {
		const Particle& part2 = Sim->particleList[(*range)[j]];

		BOOST_FOREACH(const magnet::ClonePtr<Interaction>& ptr, Sim->dynamics.getInteractions())
		  if (ptr->isInteraction(part1,part2))
		    if (dynamic_cast<const ICapture*>(ptr.get_ptr()) != NULL)
		      if (dynamic_cast<const ICapture*>(ptr.get_ptr())
			  ->isCaptured(part1,part2))
			dat.array[i * dat.chainlength + j]++;
	      }
	  }
      }
  }

  void 
  OPCContactMap::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("ContactMap");
  
    BOOST_FOREACH(Cdata& dat, chains)
      {
	//Copying
      
	for (unsigned long i = 0; i < dat.chainlength; i++)
	  {
	    for (unsigned long j = i+1; j < dat.chainlength; j++)
	      dat.array[j * dat.chainlength + i] = dat.array[i * dat.chainlength + j];
	  }

	XML << magnet::xml::tag(dat.chainPtr->getName().c_str())
	    << magnet::xml::chardata();

	//Have to draw the boxes so it renders correctly, hence the doubling up      
	for (unsigned long i = 0; i < dat.chainlength; i++)
	  {
	    for (unsigned long j = 0; j < dat.chainlength; j++)
	      XML << i - 0.5 << " " << j - 0.5 << " " << (((double) dat.array[i * dat.chainlength + j])/((double) dat.counter)) << "\n"
		  << i - 0.5 << " " << j + 0.5 << " " << (((double) dat.array[i * dat.chainlength + j])/((double) dat.counter)) << "\n";
	    XML << "\n";

	    for (unsigned long j = 0; j < dat.chainlength; j++)
	      XML << i + 0.5 << " " << j - 0.5 << " " << (((double) dat.array[i * dat.chainlength + j])/((double) dat.counter)) << "\n"
		  << i + 0.5 << " " << j + 0.5 << " " << (((double) dat.array[i * dat.chainlength + j])/((double) dat.counter)) << "\n";
	    XML << "\n";

	  }
      
            
	XML << magnet::xml::endtag(dat.chainPtr->getName().c_str());
      }

    XML << magnet::xml::endtag("ContactMap");
  }
}

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

#include "msdcorrelator.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../0partproperty/msd.hpp"
#include "../../dynamics/systems/sysTicker.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

OPMSDCorrelator::OPMSDCorrelator(const dynamo::SimData* tmp, 
				 const magnet::xml::Node& XML):
  OPTicker(tmp,"MSDCorrelator"),
  length(20),
  currCorrLength(0),
  ticksTaken(0),
  notReady(true)
{
  operator<<(XML);
}

void 
OPMSDCorrelator::operator<<(const magnet::xml::Node& XML)
{
  try 
    {
      if (XML.hasAttribute("Length"))
	length = XML.getAttribute("Length").as<size_t>();
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPMSDCorrelator";
    }    

}

void 
OPMSDCorrelator::initialise()
{
  dout << "The length of the MSD correlator is " << length << std::endl;

  posHistory.resize(Sim->N, boost::circular_buffer<Vector>(length));

  currCorrLength=1;

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    posHistory[part.getID()].push_front(part.getPosition());

  speciesData.resize(Sim->dynamics.getSpecies().size(), 
		     std::vector<double>(length, 0.0));

  structData.resize(Sim->dynamics.getTopology().size(), 
		    std::vector<double>(length, 0.0));
}

void 
OPMSDCorrelator::ticker()
{
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    posHistory[part.getID()].push_front(part.getPosition());
  
  if (notReady)
    {
      if (++currCorrLength != length)
	return;
      
      notReady = false;
    }
  
  accPass();
}

void
OPMSDCorrelator::accPass()
{
  ++ticksTaken;
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *sp->getRange())
    for (size_t step(1); step < length; ++step)
      speciesData[sp->getID()][step]
	+= (posHistory[ID][step] - posHistory[ID][0]).nrm2();
  
  BOOST_FOREACH(const magnet::ClonePtr<Topology>& topo, Sim->dynamics.getTopology())
    BOOST_FOREACH(const magnet::ClonePtr<CRange>& range, topo->getMolecules())
    {
      Vector  molCOM(0,0,0);
      double molMass(0);

      BOOST_FOREACH(const size_t& ID, *range)
	{
	  double mass = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
	  molCOM += posHistory[ID][0] * mass;
	  molMass += mass;
	}

      molCOM /= molMass;

      for (size_t step(1); step < length; ++step)
	{
	  Vector  molCOM2(0,0,0);
	  
	  BOOST_FOREACH(const size_t& ID, *range)
	    molCOM2 += posHistory[ID][step] 
	    * Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
	  
	  molCOM2 /= molMass;
	  
	  structData[topo->getID()][step] += (molCOM2 - molCOM).nrm2();
	}
    }
}

void
OPMSDCorrelator::output(magnet::xml::XmlStream &XML)
{
  XML << magnet::xml::tag("MSDCorrelator")
      << magnet::xml::tag("Particles");
  
  double dt = dynamic_cast<const CSTicker&>
    (*Sim->dynamics.getSystem("SystemTicker")).getPeriod()
    / Sim->dynamics.units().unitTime();
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    {
      XML << magnet::xml::tag("Species")
	  << magnet::xml::attr("Name")
	  << sp->getName()
	  << magnet::xml::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << speciesData[sp->getID()][step] 
	  / (static_cast<double>(ticksTaken) 
	     * static_cast<double>(sp->getCount())
	     * Sim->dynamics.units().unitArea())
	    << "\n";
      
      XML << magnet::xml::endtag("Species");
    }
  
  XML << magnet::xml::endtag("Particles")
      << magnet::xml::tag("Topology");
  
  BOOST_FOREACH(const magnet::ClonePtr<Topology>& topo, 
		Sim->dynamics.getTopology())
    {
      XML << magnet::xml::tag("Structure")
	  << magnet::xml::attr("Name")
	  << topo->getName()
	  << magnet::xml::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << structData[topo->getID()][step]
	  / (static_cast<double>(ticksTaken) 
	     * static_cast<double>(topo->getMolecules().size())
	     * Sim->dynamics.units().unitArea())
	    << "\n";
	
      XML << magnet::xml::endtag("Structure");
    }
  
  XML << magnet::xml::endtag("Topology")
      << magnet::xml::endtag("MSDCorrelator");
}

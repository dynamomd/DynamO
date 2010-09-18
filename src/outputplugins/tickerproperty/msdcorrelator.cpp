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

#include "msdcorrelator.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../0partproperty/msd.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../../dynamics/systems/sysTicker.hpp"

OPMSDCorrelator::OPMSDCorrelator(const DYNAMO::SimData* tmp, 
				   const XMLNode& XML):
  OPTicker(tmp,"MSDCorrelator"),
  length(20),
  currCorrLength(0),
  ticksTaken(0),
  notReady(true)
{
  operator<<(XML);
}

void 
OPMSDCorrelator::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	length = boost::lexical_cast<size_t>
	  (XML.getAttribute("Length"));
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPMSDCorrelator";
    }    

}

void 
OPMSDCorrelator::initialise()
{
  I_cout() << "The length of the MSD correlator is " << length;

  posHistory.resize(Sim->N, boost::circular_buffer<Vector>(length));

  currCorrLength=1;

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    posHistory[part.getID()].push_front(part.getPosition());

  speciesData.resize(Sim->dynamics.getSpecies().size(), 
		     std::vector<Iflt>(length, 0.0));

  structData.resize(Sim->dynamics.getTopology().size(), 
		    std::vector<Iflt>(length, 0.0));
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
  
  BOOST_FOREACH(const ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *sp->getRange())
    for (size_t step(1); step < length; ++step)
      speciesData[sp->getID()][step]
	+= (posHistory[ID][step] - posHistory[ID][0]).nrm2();
  
  BOOST_FOREACH(const ClonePtr<Topology>& topo, Sim->dynamics.getTopology())
    BOOST_FOREACH(const ClonePtr<CRange>& range, topo->getMolecules())
    {
      Vector  molCOM(0,0,0);
      Iflt molMass(0);

      BOOST_FOREACH(const size_t& ID, *range)
	{
	  molCOM += posHistory[ID][0] 
	    * Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass();

	  molMass += Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass();
	}

      molCOM /= molMass;

      for (size_t step(1); step < length; ++step)
	{
	  Vector  molCOM2(0,0,0);
	  
	  BOOST_FOREACH(const size_t& ID, *range)
	    molCOM2 += posHistory[ID][step] 
	    * Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass();
	  
	  molCOM2 /= molMass;
	  
	  structData[topo->getID()][step] += (molCOM2 - molCOM).nrm2();
	}
    }
}

void
OPMSDCorrelator::output(xml::XmlStream &XML)
{
  XML << xml::tag("MSDCorrelator")
      << xml::tag("Particles");
  
  Iflt dt = dynamic_cast<const CSTicker&>
    (*Sim->dynamics.getSystem("SystemTicker")).getPeriod()
    / Sim->dynamics.units().unitTime();
  
  BOOST_FOREACH(const ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    {
      XML << xml::tag("Species")
	  << xml::attr("Name")
	  << sp->getName()
	  << xml::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << speciesData[sp->getID()][step] 
	  / (static_cast<Iflt>(ticksTaken) 
	     * static_cast<Iflt>(sp->getCount())
	     * Sim->dynamics.units().unitArea())
	    << "\n";
      
      XML << xml::endtag("Species");
    }
  
  XML << xml::endtag("Particles")
      << xml::tag("Topology");
  
  BOOST_FOREACH(const ClonePtr<Topology>& topo, 
		Sim->dynamics.getTopology())
    {
      XML << xml::tag("Structure")
	  << xml::attr("Name")
	  << topo->getName()
	  << xml::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << structData[topo->getID()][step]
	  / (static_cast<Iflt>(ticksTaken) 
	     * static_cast<Iflt>(topo->getMolecules().size())
	     * Sim->dynamics.units().unitArea())
	    << "\n";
	
      XML << xml::endtag("Structure");
    }
  
  XML << xml::endtag("Topology")
      << xml::endtag("MSDCorrelator");
}

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

#include "msdcorrelator.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../0partproperty/msd.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../../dynamics/systems/sysTicker.hpp"

COPMSDCorrelator::COPMSDCorrelator(const DYNAMO::SimData* tmp, 
				   const XMLNode& XML):
  COPTicker(tmp,"MSDCorrelator"),
  length(20),
  currCorrLength(0),
  ticksTaken(0),
  notReady(true)
{
  operator<<(XML);
}

void 
COPMSDCorrelator::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	length = boost::lexical_cast<size_t>
	  (XML.getAttribute("Length"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPMSDCorrelator";
    }    
}

void 
COPMSDCorrelator::initialise()
{
  posHistory.resize(Sim->lN, boost::circular_buffer<CVector<> >(length));

  currCorrLength=1;

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    posHistory[part.getID()].push_front(part.getPosition());

  speciesData.resize(Sim->Dynamics.getSpecies().size(), 
		     std::vector<Iflt>(length, 0.0));

  structData.resize(Sim->Dynamics.getTopology().size(), 
		    std::vector<Iflt>(length, 0.0));
}

void 
COPMSDCorrelator::ticker()
{
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
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
COPMSDCorrelator::accPass()
{
  ++ticksTaken;
  
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *sp.getRange())
    for (size_t step(1); step < length; ++step)
      speciesData[sp.getID()][step] 
	+= (posHistory[ID][step] - posHistory[ID][0]).square();
  
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& topo, Sim->Dynamics.getTopology())
    BOOST_FOREACH(const smrtPlugPtr<CRange>& range, topo->getMolecules())
    {
      CVector<> molCOM(0);
      Iflt molMass(0);

      BOOST_FOREACH(const size_t& ID, *range)
	{
	  molCOM += posHistory[ID][0] 
	    * Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

	  molMass += Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();
	}

      molCOM /= molMass;

      for (size_t step(1); step < length; ++step)
	{
	  CVector<> molCOM2(0);
	  
	  BOOST_FOREACH(const size_t& ID, *range)
	    molCOM2 += posHistory[ID][step] 
	    * Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();
	  
	  molCOM2 /= molMass;
	  
	  structData[topo->getID()][step] += (molCOM2 - molCOM).square();
	}
    }
}

void
COPMSDCorrelator::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("MSDCorrelator") 
      << xmlw::tag("Particles");
  
  Iflt dt = dynamic_cast<const CSTicker&>
    (*Sim->Dynamics.getSystem("SystemTicker")).getPeriod()
    / Sim->Dynamics.units().unitTime();
  
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << sp.getName()
	  << xmlw::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << speciesData[sp.getID()][step] 
	  / (ticksTaken * sp.getCount() * Sim->Dynamics.units().unitArea())
	    << "\n";
      
      XML << xmlw::endtag("Species");
    }
  
  XML << xmlw::endtag("Particles")
      << xmlw::tag("Topology");
  
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& topo, 
		Sim->Dynamics.getTopology())
    {
      XML << xmlw::tag("Structure")
	  << xmlw::attr("Name")
	  << topo->getName()
	  << xmlw::chardata();
      
      for (size_t step(0); step < length; ++step)
	XML << dt * step << " "
	    << structData[topo->getID()][step]
	  / (ticksTaken * topo->getMolecules().size() 
	     * Sim->Dynamics.units().unitArea())
	    << "\n";
	
      XML << xmlw::endtag("Structure");
    }
  
  XML << xmlw::endtag("Topology")
      << xmlw::endtag("MSDCorrelator");
}

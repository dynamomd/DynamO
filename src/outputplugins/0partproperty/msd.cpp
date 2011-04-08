/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "msd.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

OPMSD::OPMSD(const DYNAMO::SimData* tmp, const magnet::xml::Node&):
  OutputPlugin(tmp,"MSD")
{}

OPMSD::~OPMSD()
{}

void
OPMSD::initialise()
{
  initPos.clear();
  initPos.resize(Sim->N);
  
  for (size_t ID = 0; ID < Sim->N; ++ID)
    initPos[ID] = Sim->particleList[ID].getPosition();
}

void
OPMSD::output(xml::XmlStream &XML)
{
  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();
  
  XML << xml::tag("MSD");
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    {
      double MSD(calcMSD(*(sp->getRange())));
      
      XML << xml::tag("Species")
	  << xml::attr("Name") << sp->getName()
	  << xml::attr("val") << MSD
	  << xml::attr("diffusionCoeff") 
	  << MSD * Sim->dynamics.units().unitTime() / Sim->dSysTime
	  << xml::endtag("Species");
    }

  if (!Sim->dynamics.getTopology().empty())
    {
      XML << xml::tag("Structures");

      BOOST_FOREACH(const magnet::ClonePtr<Topology>& topo, Sim->dynamics.getTopology())
	{
	  double MSD(calcStructMSD(*topo));

	  XML << xml::tag("Structure")
	      << xml::attr("Name") << topo->getName()
	      << xml::attr("val") << MSD
	      << xml::attr("diffusionCoeff") 
	      << MSD * Sim->dynamics.units().unitTime() / Sim->dSysTime
	      << xml::endtag("Structure");
	}

      XML << xml::endtag("Structures");
    }

  XML << xml::endtag("MSD");
}

double
OPMSD::calcMSD(const CRange& range) const
{
  double acc = 0.0;

  BOOST_FOREACH(const size_t ID, range)
    acc += (Sim->particleList[ID].getPosition() - initPos[ID]).nrm2();
  
  return acc / (range.size() * 2.0 * NDIM * Sim->dynamics.units().unitArea());
}

double
OPMSD::calcStructMSD(const Topology& Itop) const
{
  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();

  double acc = 0.0;
  BOOST_FOREACH(const magnet::ClonePtr<CRange>& molRange, Itop.getMolecules())
    {
      Vector  origPos(0,0,0), currPos(0,0,0);
      double totmass = 0.0;
      BOOST_FOREACH(const unsigned long& ID, *molRange)
	{
	  double pmass = Sim->dynamics.getSpecies(Sim->particleList[ID])
	    .getMass();

	  totmass += pmass;
	  currPos += Sim->particleList[ID].getPosition() * pmass;
	  origPos += initPos[ID] * pmass;
	}
      
      currPos /= totmass;
      origPos /= totmass;

      acc += (currPos - origPos).nrm2();
    }

  acc /= Itop.getMoleculeCount() * 2.0 * NDIM
    * Sim->dynamics.units().unitArea();
  
  return acc;
}

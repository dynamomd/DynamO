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

#include "msd.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"


OPMSD::OPMSD(const DYNAMO::SimData* tmp, const XMLNode&):
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
  {
    Iflt MSD(calcMSD());

    XML << xml::tag("MSD") 
	<< xml::tag("Particle") 
	<< xml::attr("val") << MSD
	<< xml::attr("diffusionCoeff") 
	<< MSD * Sim->dynamics.units().unitTime() / Sim->dSysTime
	<< xml::endtag("Particle");
  }

  if (!Sim->dynamics.getTopology().empty())
    {
      XML << xml::tag("Structures");

      BOOST_FOREACH(const ClonePtr<Topology>& topo, Sim->dynamics.getTopology())
	{
	  Iflt MSD(calcStructMSD(*topo));

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

Iflt
OPMSD::calcMSD() const
{
  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();

  Iflt acc = 0.0;
  
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    acc += (part.getPosition() - initPos[part.getID()]).nrm2();
  
  return acc / (initPos.size() * 2.0 * NDIM * Sim->dynamics.units().unitArea());
}

Iflt
OPMSD::calcStructMSD(const Topology& Itop) const
{
  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();

  Iflt acc = 0.0;
  BOOST_FOREACH(const ClonePtr<CRange>& molRange, Itop.getMolecules())
    {
      Vector  origPos(0,0,0), currPos(0,0,0);
      Iflt totmass = 0.0;
      BOOST_FOREACH(const unsigned long& ID, *molRange)
	{
	  Iflt pmass = Sim->dynamics.getSpecies(Sim->particleList[ID])
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

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

#include "msd.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"


COPMSD::COPMSD(const DYNAMO::SimData* tmp, const XMLNode&):
  COutputPlugin(tmp,"MSD")
{}

COPMSD::~COPMSD()
{}

void
COPMSD::initialise()
{
  initPos.clear();
  initPos.resize(Sim->lN);
  
  for (size_t ID = 0; ID < Sim->lN; ++ID)
    initPos[ID] = Sim->vParticleList[ID].getPosition();
}

void
COPMSD::output(xmlw::XmlStream &XML)
{
  {
    Iflt MSD(calcMSD());

    XML << xmlw::tag("MSD") 
	<< xmlw::tag("Particle") 
	<< xmlw::attr("val") << MSD
	<< xmlw::attr("diffusionCoeff") 
	<< MSD * Sim->Dynamics.units().unitTime() / Sim->dSysTime
	<< xmlw::endtag("Particle");
  }

  if (!Sim->Dynamics.getTopology().empty())
    {
      XML << xmlw::tag("Structures");

      BOOST_FOREACH(const smrtPlugPtr<CTopology>& topo, Sim->Dynamics.getTopology())
	{
	  Iflt MSD(calcStructMSD(*topo));

	  XML << xmlw::tag("Structure")
	      << xmlw::attr("Name") << topo->getName()
	      << xmlw::attr("val") << MSD
	      << xmlw::attr("diffusionCoeff") 
	      << MSD * Sim->Dynamics.units().unitTime() / Sim->dSysTime
	      << xmlw::endtag("Structure");
	}

      XML << xmlw::endtag("Structures");
    }

  XML << xmlw::endtag("MSD");
}

Iflt
COPMSD::calcMSD() const
{
  //Required to get the correct results
  Sim->Dynamics.Liouvillean().updateAllParticles();

  Iflt acc = 0.0;
  
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    acc += (part.getPosition() - initPos[part.getID()]).square();
  
  return acc / (initPos.size() * 2.0 * NDIM * Sim->Dynamics.units().unitArea());
}

Iflt
COPMSD::calcStructMSD(const CTopology& Itop) const
{
  //Required to get the correct results
  Sim->Dynamics.Liouvillean().updateAllParticles();

  Iflt acc = 0.0;
  BOOST_FOREACH(const smrtPlugPtr<CRange>& molRange, Itop.getMolecules())
    {
      CVector<> origPos(0), currPos(0);
      Iflt totmass = 0.0;
      BOOST_FOREACH(const unsigned long& ID, *molRange)
	{
	  Iflt pmass = Sim->Dynamics.getSpecies(Sim->vParticleList[ID])
	    .getMass();

	  totmass += pmass;
	  currPos += Sim->vParticleList[ID].getPosition() * pmass;
	  origPos += initPos[ID] * pmass;
	}
      
      currPos /= totmass;
      origPos /= totmass;

      acc += (currPos - origPos).square();
    }

  acc /= Itop.getMoleculeCount() * 2.0 * NDIM
    * Sim->Dynamics.units().unitArea();
  
  return acc;
}

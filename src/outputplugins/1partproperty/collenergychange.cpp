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

#include "collenergychange.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_simdata.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/species/species.hpp"
#include "../../dynamics/1particleEventData.hpp"
#include "../../dynamics/2particleEventData.hpp"
#include "../../dynamics/units/units.hpp"

OPCollEnergyChange::OPCollEnergyChange(const DYNAMO::SimData* tmp, const XMLNode&XML):
  OP1PP(tmp,"CollEnergyChange", 250),
  binWidth(0.001)
{ operator<<(XML); }

void 
OPCollEnergyChange::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("binWidth"))
      binWidth = boost::lexical_cast<Iflt>(XML.getAttribute("binWidth"));
      }
  catch (std::exception& excep)
    {
      D_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void
OPCollEnergyChange::initialise()
{
  I_cout() << "Bin width set to " << binWidth;

  data.resize(Sim->dynamics.getSpecies().size(), 
	      C1DHistogram(Sim->dynamics.units().unitEnergy() * binWidth));

  specialhist = C1DHistogram(Sim->dynamics.units().unitEnergy() * binWidth);
}

void 
OPCollEnergyChange::A1ParticleChange(const C1ParticleData& PDat)
{
  data[PDat.getSpecies().getID()]
    .addVal(PDat.getDeltaKE());
}

void 
OPCollEnergyChange::A2ParticleChange(const C2ParticleData& PDat)
{
  data[PDat.particle1_.getSpecies().getID()]
    .addVal(PDat.particle1_.getDeltaKE());

  data[PDat.particle2_.getSpecies().getID()]
    .addVal(PDat.particle2_.getDeltaKE());

  Iflt p1Mass = PDat.particle1_.getSpecies().getMass(); 
  Iflt p2Mass = PDat.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass / (p1Mass+p2Mass);

  specialhist.addVal((PDat.dP.nrm2() / (2.0 * mu)) - (PDat.vijold | PDat.dP));
}

void
OPCollEnergyChange::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("CollEnergyChange")
      << xmlw::tag("PairCalc");

  specialhist.outputHistogram(XML, 1.0 / Sim->dynamics.units().unitEnergy());

  XML << xmlw::endtag("PairCalc");

  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->dynamics.getSpecies()[id]->getName();

      data[id].outputHistogram(XML, 1.0 / Sim->dynamics.units().unitEnergy());

      XML << xmlw::endtag("Species");
    }

  XML << xmlw::endtag("CollEnergyChange");
}


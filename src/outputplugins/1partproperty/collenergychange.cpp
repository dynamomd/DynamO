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

#include "collenergychange.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/species/species.hpp"
#include "../../dynamics/1particleEventData.hpp"
#include "../../dynamics/2particleEventData.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../extcode/xmlwriter.hpp"

double OPCollEnergyChange::KEBinWidth = 0.01;


OPCollEnergyChange::OPCollEnergyChange(const DYNAMO::SimData* tmp, const XMLNode&XML):
  OP1PP(tmp,"CollEnergyChange", 250),
  binWidth(0.001)
{ operator<<(XML); }

void 
OPCollEnergyChange::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("binWidth"))
      binWidth = boost::lexical_cast<double>(XML.getAttribute("binWidth"));

    if (XML.isAttributeSet("KEBinWidth"))
      KEBinWidth = boost::lexical_cast<double>(XML.getAttribute("KEBinWidth"));

    KEBinWidth *= Sim->dynamics.units().unitEnergy();
      }
  catch (std::exception& excep)
    {
      M_throw() << "Error while parsing " << name << "options\n"
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
OPCollEnergyChange::A1ParticleChange(const ParticleEventData& PDat)
{
  data[PDat.getSpecies().getID()]
    .addVal(PDat.getDeltaKE());
}

void 
OPCollEnergyChange::A2ParticleChange(const PairEventData& PDat)
{
  data[PDat.particle1_.getSpecies().getID()]
    .addVal(PDat.particle1_.getDeltaKE());

  data[PDat.particle2_.getSpecies().getID()]
    .addVal(PDat.particle2_.getDeltaKE());

  double p1Mass = PDat.particle1_.getSpecies().getMass(); 
  double p2Mass = PDat.particle2_.getSpecies().getMass();
  double mu = p1Mass * p2Mass / (p1Mass+p2Mass);

  specialhist.addVal((PDat.dP.nrm2() / (2.0 * mu)) - (PDat.vijold | PDat.dP));

  collisionKE[mapkey(PDat.particle1_.getSpecies().getID(), 
		     PDat.particle2_.getSpecies().getID(), 
		     PDat.getType())]
    .addVal(Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.particle1_.getParticle())
	    -PDat.particle1_.getDeltaKE());

  collisionKE[mapkey(PDat.particle2_.getSpecies().getID(), 
		     PDat.particle1_.getSpecies().getID(), 
		     PDat.getType())]
    .addVal(Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.particle2_.getParticle())
	    -PDat.particle2_.getDeltaKE());
}

void
OPCollEnergyChange::output(xml::XmlStream &XML)
{
  XML << xml::tag("CollEnergyChange")
      << xml::tag("PairCalc");

  specialhist.outputHistogram(XML, 1.0 / Sim->dynamics.units().unitEnergy());

  XML << xml::endtag("PairCalc");

  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xml::tag("Species")
	  << xml::attr("Name")
	  << Sim->dynamics.getSpecies()[id]->getName();

      data[id].outputHistogram(XML, 1.0 / Sim->dynamics.units().unitEnergy());

      XML << xml::endtag("Species");
    }

  typedef std::pair<const mapkey, histogram> locpair;
  BOOST_FOREACH(const locpair& pdat, collisionKE)
    {
      XML << xml::tag("Energy_On_Collision")
	  << xml::attr("Species") << Sim->dynamics.getSpecies()[pdat.first.get<0>()]->getName()
	  << xml::attr("EventPartnerSpecies")
	  << Sim->dynamics.getSpecies()[pdat.first.get<1>()]->getName()
	  << xml::attr("EventType") 
	  << pdat.first.get<2>();
      
      pdat.second.outputHistogram(XML, 1.0 / Sim->dynamics.units().unitEnergy());      

      XML << xml::endtag("Energy_On_Collision");
    }

  XML << xml::endtag("CollEnergyChange");
}


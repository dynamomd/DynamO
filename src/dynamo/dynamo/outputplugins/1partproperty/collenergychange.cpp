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

#include <dynamo/outputplugins/1partproperty/collenergychange.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>

#include <dynamo/species/species.hpp>
#include <dynamo/1particleEventData.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  double OPCollEnergyChange::KEBinWidth = 0.01;


  OPCollEnergyChange::OPCollEnergyChange(const dynamo::Simulation* tmp, 
					 const magnet::xml::Node&XML):
    OP1PP(tmp,"CollEnergyChange", 250),
    binWidth(0.001)
  { operator<<(XML); }

  void 
  OPCollEnergyChange::operator<<(const magnet::xml::Node& XML)
  {
    try {
      if (XML.hasAttribute("binWidth"))
	binWidth = XML.getAttribute("binWidth").as<double>();

      if (XML.hasAttribute("KEBinWidth"))
	KEBinWidth = XML.getAttribute("KEBinWidth").as<double>();
      KEBinWidth *= Sim->units.unitEnergy();
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
    dout << "Bin width set to " << binWidth << std::endl;

    data.resize(Sim->species.size(), 
		magnet::math::Histogram<>(Sim->units.unitEnergy() * binWidth));

    specialhist = magnet::math::Histogram<>(Sim->units.unitEnergy() * binWidth);
  }

  void 
  OPCollEnergyChange::A1ParticleChange(const ParticleEventData& PDat)
  {
    data[PDat.getSpeciesID()].addVal(PDat.getDeltaKE());
  }

  void 
  OPCollEnergyChange::A2ParticleChange(const PairEventData& PDat)
  {
    data[PDat.particle1_.getSpeciesID()].addVal(PDat.particle1_.getDeltaKE());
    data[PDat.particle2_.getSpeciesID()].addVal(PDat.particle2_.getDeltaKE());

    const Particle& p1 = Sim->particles[PDat.particle1_.getParticleID()];
    const Particle& p2 = Sim->particles[PDat.particle2_.getParticleID()];

    double p1Mass = Sim->species[PDat.particle1_.getSpeciesID()]->getMass(p1.getID()); 
    double p2Mass = Sim->species[PDat.particle2_.getSpeciesID()]->getMass(p2.getID());
    double mu = p1Mass * p2Mass / (p1Mass + p2Mass);

    specialhist.addVal((PDat.impulse.nrm2() / (2.0 * mu)) - (PDat.vijold | PDat.impulse));

    collisionKE[mapkey(PDat.particle1_.getSpeciesID(), 
		       PDat.particle2_.getSpeciesID(), 
		       PDat.getType())]
      .addVal(Sim->dynamics->getParticleKineticEnergy(p1)
	      -PDat.particle1_.getDeltaKE());

    collisionKE[mapkey(PDat.particle2_.getSpeciesID(), 
		       PDat.particle1_.getSpeciesID(), 
		       PDat.getType())]
      .addVal(Sim->dynamics->getParticleKineticEnergy(p2)
	      -PDat.particle2_.getDeltaKE());
  }

  void
  OPCollEnergyChange::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("CollEnergyChange")
	<< magnet::xml::tag("PairCalc");

    specialhist.outputHistogram(XML, 1.0 / Sim->units.unitEnergy());

    XML << magnet::xml::endtag("PairCalc");

    for (size_t id = 0; id < data.size(); ++id)
      {
	XML << magnet::xml::tag("Species")
	    << magnet::xml::attr("Name")
	    << Sim->species[id]->getName();

	data[id].outputHistogram(XML, 1.0 / Sim->units.unitEnergy());

	XML << magnet::xml::endtag("Species");
      }

    typedef std::pair<const mapkey, histogram> locpair;
    BOOST_FOREACH(const locpair& pdat, collisionKE)
      {
	XML << magnet::xml::tag("Energy_On_Collision")
	    << magnet::xml::attr("Species") << Sim->species[pdat.first.get<0>()]->getName()
	    << magnet::xml::attr("EventPartnerSpecies")
	    << Sim->species[pdat.first.get<1>()]->getName()
	    << magnet::xml::attr("EventType") 
	    << pdat.first.get<2>();
      
	pdat.second.outputHistogram(XML, 1.0 / Sim->units.unitEnergy());      

	XML << magnet::xml::endtag("Energy_On_Collision");
      }

    XML << magnet::xml::endtag("CollEnergyChange");
  }
}

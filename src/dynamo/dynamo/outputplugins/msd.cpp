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

#include <dynamo/include.hpp>
#include <dynamo/outputplugins/msd.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPMSD::OPMSD(const dynamo::Simulation *tmp, const magnet::xml::Node &)
    : OutputPlugin(tmp, "MSD") {}

OPMSD::~OPMSD() {}

void OPMSD::initialise() {
  initPos.clear();
  initPos.resize(Sim->N());

  for (size_t ID = 0; ID < Sim->N(); ++ID)
    initPos[ID] = Sim->particles[ID].getPosition();
}

void OPMSD::output(magnet::xml::XmlStream &XML) {
  // Required to get the correct results
  Sim->dynamics->updateAllParticles();

  XML << magnet::xml::tag("MSD");

  for (const shared_ptr<Species> &sp : Sim->species) {
    Vector MSD = calcMSD(*sp->getRange()) / Sim->units.unitArea();
    double MSD_sum = 0;
    for (size_t i(0); i < NDIM; ++i)
      MSD_sum += MSD[i];
    MSD_sum /= 3;
    XML << magnet::xml::tag("Species") << magnet::xml::attr("Name")
        << sp->getName() << magnet::xml::attr("val") << MSD_sum
        << magnet::xml::attr("diffusionCoeff")
        << MSD_sum * Sim->units.unitTime() / (2 * Sim->systemTime)
        << magnet::xml::tag("MSDvec") << MSD << magnet::xml::endtag("MSDvec")
        << magnet::xml::tag("Dvec")
        << MSD * Sim->units.unitTime() / (2 * Sim->systemTime)
        << magnet::xml::endtag("Dvec") << magnet::xml::endtag("Species");
  }

  if (!Sim->topology.empty()) {
    XML << magnet::xml::tag("Structures");

    for (const shared_ptr<Topology> &topo : Sim->topology) {
      Vector MSD = calcStructMSD(*topo) / Sim->units.unitArea();

      double MSD_sum = 0;
      for (size_t i(0); i < NDIM; ++i)
        MSD_sum += MSD[i];
      MSD_sum /= 3;

      XML << magnet::xml::tag("Structure") << magnet::xml::attr("Name")
          << topo->getName() << magnet::xml::attr("val") << MSD_sum
          << magnet::xml::attr("diffusionCoeff")
          << MSD_sum * Sim->units.unitTime() / (2 * NDIM * Sim->systemTime)
          << magnet::xml::tag("MSDvec") << MSD << magnet::xml::endtag("MSDvec")
          << magnet::xml::tag("Dvec")
          << MSD * Sim->units.unitTime() / (2 * Sim->systemTime)
          << magnet::xml::endtag("Dvec") << magnet::xml::endtag("Structure");
    }

    XML << magnet::xml::endtag("Structures");
  }

  XML << magnet::xml::endtag("MSD");
}

Vector OPMSD::calcMSD(const IDRange &range) const {
  Vector acc(0.0);

  for (const size_t ID : range) {
    const Vector diff = Sim->particles[ID].getPosition() - initPos[ID];
    for (size_t i(0); i < NDIM; ++i)
      acc[i] += diff[i] * diff[i];
  }

  return acc / range.size();
}

Vector OPMSD::calcD(const IDRange &range) const {
  return calcMSD(range) / (2 * Sim->systemTime);
}

Vector OPMSD::calcStructMSD(const Topology &Itop) const {
  // Required to get the correct results
  Sim->dynamics->updateAllParticles();

  Vector acc(0.0);
  for (const shared_ptr<IDRange> &molRange : Itop.getMolecules()) {
    Vector origPos{0, 0, 0}, currPos{0, 0, 0};
    double totmass = 0.0;
    for (const unsigned long &ID : *molRange) {
      double pmass = Sim->species[Sim->particles[ID]]->getMass(ID);

      totmass += pmass;
      currPos += Sim->particles[ID].getPosition() * pmass;
      origPos += initPos[ID] * pmass;
    }

    currPos /= totmass;
    origPos /= totmass;
    const Vector diff = currPos - origPos;
    for (size_t i(0); i < NDIM; ++i)
      acc[i] += diff[i] * diff[i];
  }

  acc /= Itop.getMoleculeCount();

  return acc;
}
} // namespace dynamo

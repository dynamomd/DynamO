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

#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/include.hpp>
#include <dynamo/outputplugins/tickerproperty/vacf.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/sysTicker.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPVACF::OPVACF(const dynamo::Simulation *tmp, const magnet::xml::Node &XML)
    : OPTicker(tmp, "VACF"), length(50), currCorrLength(0), ticksTaken(0) {
  OPVACF::operator<<(XML);
}

void OPVACF::operator<<(const magnet::xml::Node &XML) {
  if (XML.hasAttribute("Length"))
    length = XML.getAttribute("Length").as<size_t>();
}

void OPVACF::initialise() {
  dout << "The length of the VACF correlator is " << length << std::endl;

  velHistory.resize(Sim->N(), boost::circular_buffer<Vector>(length));

  currCorrLength = 1;

  for (const Particle &part : Sim->particles)
    velHistory[part.getID()].push_front(part.getPosition());

  speciesData.resize(Sim->species.size(), std::vector<double>(length, 0.0));
  structData.resize(Sim->topology.size(), std::vector<double>(length, 0.0));
}

void OPVACF::ticker() {
  for (const Particle &part : Sim->particles)
    velHistory[part.getID()].push_front(part.getVelocity());

  currCorrLength +=
      (currCorrLength !=
       length); // Add to the length while we're not at full length

  accPass();
}

void OPVACF::accPass() {
  ++ticksTaken;

  for (const shared_ptr<Species> &sp : Sim->species)
    for (const size_t &ID : *sp->getRange())
      for (size_t step(0); step < currCorrLength; ++step)
        speciesData[sp->getID()][step] +=
            velHistory[ID][step] | velHistory[ID][0];

  for (const shared_ptr<Topology> &topo : Sim->topology)
    for (const shared_ptr<IDRange> &range : topo->getMolecules()) {
      Vector COMvelocity({0, 0, 0});
      double molMass(0);

      for (const size_t &ID : *range) {
        const auto &p = Sim->particles[ID];
        double mass = (Sim->species[p])->getMass(ID);
        COMvelocity += velHistory[ID][0] * mass;
        molMass += mass;
      }

      COMvelocity /= molMass;

      for (size_t step(0); step < currCorrLength; ++step) {
        Vector COMvelocity2({0, 0, 0});

        for (const size_t &ID : *range)
          COMvelocity2 += velHistory[ID][step] *
                          Sim->species[Sim->particles[ID]]->getMass(ID);
        COMvelocity2 /= molMass;
        structData[topo->getID()][step] += COMvelocity | COMvelocity2;
      }
    }
}

void OPVACF::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("VACF") << magnet::xml::attr("ticks") << ticksTaken
      << magnet::xml::tag("Particles");

  const double dt =
      dynamic_cast<const SysTicker &>(*Sim->systems["SystemTicker"])
          .getPeriod() /
      Sim->units.unitTime();

  for (const shared_ptr<Species> &sp : Sim->species) {
    XML << magnet::xml::tag("Species") << magnet::xml::attr("Name")
        << sp->getName() << magnet::xml::chardata();

    for (size_t step(0); step < currCorrLength; ++step)
      XML << dt * step << " "
          << speciesData[sp->getID()][step] /
                 (static_cast<double>(ticksTaken - step) *
                  static_cast<double>(sp->getCount()) *
                  Sim->units.unitVelocity() * Sim->units.unitVelocity())
          << "\n";

    XML << magnet::xml::endtag("Species");
  }

  XML << magnet::xml::endtag("Particles") << magnet::xml::tag("Topology");

  for (const shared_ptr<Topology> &topo : Sim->topology) {
    XML << magnet::xml::tag("Structure") << magnet::xml::attr("Name")
        << topo->getName() << magnet::xml::chardata();

    for (size_t step(0); step < currCorrLength; ++step)
      XML << dt * step << " "
          << structData[topo->getID()][step] /
                 (static_cast<double>(ticksTaken - step) *
                  static_cast<double>(topo->getMolecules().size()) *
                  Sim->units.unitVelocity() * Sim->units.unitVelocity())
          << "\n";

    XML << magnet::xml::endtag("Structure");
  }

  XML << magnet::xml::endtag("Topology") << magnet::xml::endtag("VACF");
}
} // namespace dynamo

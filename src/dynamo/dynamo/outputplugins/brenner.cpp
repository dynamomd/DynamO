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
#include <dynamo/outputplugins/brenner.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPBrenner::OPBrenner(const dynamo::Simulation *tmp, const magnet::xml::Node &)
    : OutputPlugin(tmp, "Brenner") {}

OPBrenner::~OPBrenner() {}

void OPBrenner::initialise() {
  for (size_t i(0); i < NDIM; ++i)
    _sysmomentum_hist[i] = magnet::math::HistogramWeighted<>(
        Sim->N() * 0.001 * Sim->units.unitMomentum());

  _sysMomentum = Vector({0, 0, 0});

  for (const Particle &part : Sim->particles) {
    const Species &sp = *(Sim->species[part]);
    const double mass = sp.getMass(part.getID());
    if (std::isinf(mass))
      continue;
    _sysMomentum += mass * part.getVelocity();
  }
}

void OPBrenner::eventUpdate(const Event &event, const NEventData &SDat) {
  for (size_t i(0); i < NDIM; ++i)
    _sysmomentum_hist[i].addVal(_sysMomentum[i], event._dt);

  for (const ParticleEventData &pData : SDat.L1partChanges) {
    const Particle &p1 = Sim->particles[pData.getParticleID()];
    const double m1 = Sim->species[p1]->getMass(p1.getID());
    const Vector dP = m1 * (p1.getVelocity() - pData.getOldVel());
    _sysMomentum += dP;
  }
}

void OPBrenner::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("OPBrenner");
  for (size_t i(0); i < NDIM; ++i) {
    XML << magnet::xml::tag("Dimension") << magnet::xml::attr("dim") << i;
    _sysmomentum_hist[i].outputClearHistogram(XML, Sim->units.unitMomentum());
    XML << magnet::xml::endtag("Dimension");
  }

  XML << magnet::xml::endtag("OPBrenner");
}
} // namespace dynamo

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
#include <dynamo/outputplugins/tickerproperty/velprof.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPVelProfile::OPVelProfile(const dynamo::Simulation *tmp,
                           const magnet::xml::Node &)
    : OPTicker(tmp, "VelProfile"), samplesTaken(0), binWidth(0.5) {
  if (NDIM != 3)
    M_throw() << "Terrible plugin for 3 dims only";
}

void OPVelProfile::initialise() {
  binWidth *= Sim->units.unitLength();

  vx.resize(static_cast<size_t>(Sim->primaryCellSize[1] / binWidth) + 1,
            std::vector<std::pair<size_t, double>>(
                static_cast<size_t>(Sim->primaryCellSize[2] / binWidth) + 1,
                std::pair<size_t, double>(0, 0)));
}

void OPVelProfile::ticker() {
  for (const Particle &part : Sim->particles) {
    Vector pos(part.getPosition());
    Vector vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    pos += Sim->primaryCellSize / 2.0;

    std::pair<size_t, double> &locpair(vx[static_cast<size_t>(
        pos[1] / binWidth)][static_cast<size_t>(pos[2] / binWidth)]);
    ++locpair.first;
    locpair.second += vel[0];
  }
}

void OPVelProfile::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("VelProfile") << magnet::xml::chardata();

  size_t nybins = static_cast<size_t>(Sim->primaryCellSize[1] / binWidth) + 1;
  size_t nzbins = static_cast<size_t>(Sim->primaryCellSize[2] / binWidth) + 1;

  for (size_t y = 0; y < nybins; ++y) {
    for (size_t z = 0; z < nzbins; ++z)
      XML << y * binWidth / Sim->units.unitLength() << " "
          << z * binWidth / Sim->units.unitLength() << " "
          << ((vx[y][z].first) ? (vx[y][z].second / vx[y][z].first) /
                                     Sim->units.unitVelocity()
                               : 0)
          << "\n";

    XML << "\n";
  }

  XML << magnet::xml::endtag("VelProfile");
}
} // namespace dynamo

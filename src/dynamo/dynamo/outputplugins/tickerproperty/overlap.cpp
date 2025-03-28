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
#include <dynamo/outputplugins/tickerproperty/overlap.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
OPOverlapTest::OPOverlapTest(const dynamo::Simulation *tmp,
                             const magnet::xml::Node &)
    : OPTicker(tmp, "OverlapTester") {}

void OPOverlapTest::initialise() {
  dout << "Testing for overlaps in starting configuration" << std::endl;
  ticker();
}

void OPOverlapTest::output(magnet::xml::XmlStream &) {
  dout << "Testing for overlaps in output configuration" << std::endl;
  ticker();
}

void OPOverlapTest::ticker() {
  for (std::vector<Particle>::const_iterator iPtr = Sim->particles.begin();
       iPtr != Sim->particles.end(); ++iPtr)
    for (std::vector<Particle>::const_iterator jPtr = iPtr + 1;
         jPtr != Sim->particles.end(); ++jPtr)
      Sim->getInteraction(*iPtr, *jPtr)->validateState(*iPtr, *jPtr);
}
} // namespace dynamo

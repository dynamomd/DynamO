
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
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/systems/sysTicker.hpp>

namespace dynamo {
OPTicker::OPTicker(const dynamo::Simulation *t1, const char *t2)
    : OutputPlugin(t1, t2) {}

double OPTicker::getTickerTime() const {
  try {
    return std::dynamic_pointer_cast<SysTicker>(Sim->systems["SystemTicker"])
        ->getPeriod();
  } catch (const std::bad_cast &) {
    M_throw() << "Could not upcast the SystemTicker system event to SysTicker, "
                 "have you named a system as SystemTicker?";
  }
}
} // namespace dynamo

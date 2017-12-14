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

#include <dynamo/base.hpp>
#include <magnet/stream/console_specials.hpp>

namespace dynamo
{
  Base::Base(const std::string name):
    dout(std::cout), derr(std::cerr)
  {
    setOutputPrefix(name);
    //Reasonable precision for output
    dout << std::setprecision(std::numeric_limits<float>::digits10);
    derr << std::setprecision(std::numeric_limits<float>::digits10);
  }

  void Base::setOutputPrefix(const std::string& prefix)
  {
#ifdef DYNAMO_COLORIZE
    dout.setPrefix(colorCode(prefix) + prefix + ": " + magnet::console::reset());
    derr.setPrefix(magnet::console::bold() + magnet::console::red_fg() + prefix + ": " + magnet::console::reset());
#else
    dout.setPrefix(prefix + ": ");
    derr.setPrefix(prefix + ": ");
#endif
  }

  Base::Base() { M_throw() << "Calling the default constructor!"; }

  std::string
  Base::colorCode(std::string str)
  {
    switch (std::hash<std::string>()(str) % 9)
      {
      case 0: return magnet::console::cyan_fg();
      case 1: return magnet::console::purple_fg();
      case 2: return magnet::console::blue_fg();
      case 3: return magnet::console::yellow_fg();
      case 4: return magnet::console::green_fg();
      case 5: return magnet::console::bold() + magnet::console::green_fg();
      case 6: return magnet::console::bold() + magnet::console::blue_fg();
      case 7: return magnet::console::bold() + magnet::console::purple_fg();
      case 8: return magnet::console::bold() + magnet::console::cyan_fg();
      }
    return "";
  }

  SimBase::SimBase(Simulation* const SD, const std::string aName):
    Base(aName),
    Sim(SD)
  {}

  SimBase_const::SimBase_const(const Simulation* const SD, const std::string aName):
    Base(aName),
    Sim(SD)
  {}
}

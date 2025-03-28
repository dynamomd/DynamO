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

#pragma once
#include <array>
#include <boost/program_options.hpp>
#include <dynamo/base.hpp>
#include <magnet/math/vector.hpp>

using namespace std;
using namespace boost;
namespace po = boost::program_options;

namespace dynamo {
class UCell;

class IPPacker : public dynamo::SimBase {
public:
  IPPacker(po::variables_map &, dynamo::Simulation *tmp);

  void initialise();

  static po::options_description getOptions();

protected:
  std::array<long, 3> getCells();
  Vector getNormalisedCellDimensions();
  Vector getRandVelVec();
  UCell *standardPackingHelper(UCell *, bool forceRectangular = false);

  po::variables_map &vm;
};
} // namespace dynamo

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
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/systems/sysTicker.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace dynamo {

  Dynamics::Dynamics(dynamo::SimData* tmp): 
    SimBase(tmp, "Dynamics")
  {}
}

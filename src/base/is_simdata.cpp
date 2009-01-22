/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "is_simdata.hpp"
#include "../schedulers/scheduler.hpp"

namespace DYNAMO
{
  SimData::SimData():
    Ensemble(NULL),
    dSysTime(0.0),
    lNColl(0),
    lMaxNColl(100000),
    lNPrint(50000),
    lPrintLimiter(0),
    lN(0),
    ptrScheduler(NULL),
    Dynamics(this),
    aspectRatio(1.0),
    ranGenerator(static_cast<unsigned>(std::time(0))),
    normal_sampler(ranGenerator, boost::normal_distribution_01<Iflt>()),
    uniform_sampler(ranGenerator),
    lastRunMFT(0.0),
    simID(0),
    status(START),
    binaryXML(false)
  {}

  SimData::~SimData()
  {
    if (ptrScheduler != NULL) delete ptrScheduler;
  }
}

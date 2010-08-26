/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "eta.hpp"
#include "../../base/is_simdata.hpp"

OPETA::OPETA(const DYNAMO::SimData* tmp, const XMLNode&):
  OutputPlugin(tmp,"EstTime", 249)
{}

void
OPETA::initialise()
{
  time(&start_Time);
}

void
OPETA::periodicOutput()
{
  time_t currTime;
  time(&currTime);

  double ETA = floor(((Sim->endEventCount - Sim->eventCount) * difftime(currTime, start_Time)) / Sim->eventCount);
  double ETA_hours = floor(ETA/3600);
  double ETA_mins = floor(ETA/60) - (ETA_hours * 60);
  double ETA_secs = ETA - (ETA_mins * 60) - (ETA_hours * 3600);

  I_Pcout() << "ETA " << ETA_hours << "h " << ETA_mins << "m " << ETA_secs << "s, ";
}

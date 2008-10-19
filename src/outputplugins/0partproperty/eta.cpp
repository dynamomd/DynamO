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

#include "eta.hpp"
#include "../../base/is_simdata.hpp"

COPETA::COPETA(const DYNAMO::SimData* tmp, const XMLNode&):
  COutputPlugin(tmp,"EstTime", 249)
{}

void
COPETA::initialise() 
{
  time(&start_Time);
}

void
COPETA::periodicOutput()
{
  time_t currTime;
  time(&currTime);
  
  I_Pcout() << "ETA " << ((Sim->lMaxNColl - Sim->lNColl) 
			  * difftime(currTime, start_Time)) / Sim->lNColl  
	    << "s, ";
}

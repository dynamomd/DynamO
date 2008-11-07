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

#include "Null_sched.hpp"
#include "../base/is_exception.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/globEvent.hpp"

CSNull::CSNull(const DYNAMO::SimData* Sim):
  CScheduler(Sim,"NullSched")
{
  I_cout() << "Null scheduler loaded";
}

void 
CSNull::initialise() 
{}

void 
CSNull::update(const CParticle&) 
{}

void 
CSNull::stream(const Iflt&) {}

ENextEvent 
CSNull::nextEventType() const
{ D_throw() << "CSNull is not a real scheduler!"; }

void 
CSNull::rebuildList() 
{}


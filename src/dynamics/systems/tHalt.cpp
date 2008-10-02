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

#include "tHalt.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../units/units.hpp"

CStHalt::CStHalt(DYNAMO::SimData* nSim, Iflt ndt, std::string nName):
  CSystem(nSim)
{
  dt = ndt * Sim->Dynamics.units().unitTime();

  sysName = nName;

  I_cout() << "System halt set for " 
	   << ndt;
}

void 
CStHalt::stream(Iflt ndt)
{
  dt -= ndt;
}

CNParticleData 
CStHalt::runEvent()
{
  Sim->lPrintLimiter = Sim->lMaxNColl = Sim->lNColl;
  return CNParticleData();
}

void 
CStHalt::initialise(size_t nID)
{
  ID=nID;
}

void 
CStHalt::setdt(Iflt ndt)
{ dt = ndt * Sim->Dynamics.units().unitTime(); }

void 
CStHalt::increasedt(Iflt ndt)
{ 
  dt += ndt * Sim->Dynamics.units().unitTime(); 
}

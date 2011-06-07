/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "nblistCompressionFix.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/neighbourList.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../../schedulers/scheduler.hpp"

CSNBListCompressionFix::CSNBListCompressionFix(dynamo::SimData* nSim, double nGR, size_t nblistID):
  System(nSim),
  growthRate(nGR),
  cellID(nblistID)
{
  sysName = "GlobalCellsCompressionHack";
  type = NON_EVENT;

  if (dynamic_cast<const CGNeighbourList*>(Sim->dynamics.getGlobals()[cellID].get_ptr()) == NULL)
    M_throw() << "The ID passed to CSNBListCompressionFix isn't a CGNeighbourList";
}

void
CSNBListCompressionFix::initialise(size_t nID)
{
  ID = nID;

  if (dynamic_cast<const CGNeighbourList*>(Sim->dynamics.getGlobals()[cellID].get_ptr()) == NULL)
    M_throw() << "Have the globals been shuffled? The cellID is no longer a CGNeighbourList.";
  
  CGNeighbourList& nblist(dynamic_cast<CGNeighbourList&>(*Sim->dynamics.getGlobals()[cellID]));

  dt = (nblist.getMaxSupportedInteractionLength() / nblist.getMaxInteractionLength() - 1.0) / growthRate - Sim->dSysTime;

  dout << "Compression Hack Loaded"
	   << "\nFor global " << nblist.getName()
	   << "\nCompression rate = " 
	   << growthRate / Sim->dynamics.units().unitTime()
	   << "\nSim Units compression rate = " << growthRate
	   << "\nMax length of interaction = " 
	   << nblist.getMaxSupportedInteractionLength() / Sim->dynamics.units().unitLength()
	   << "\nMaximum supported length = "
	   << nblist.getMaxSupportedInteractionLength() / Sim->dynamics.units().unitLength()
	   << "\nFirst halt scheduled for " 
	   << dt / Sim->dynamics.units().unitTime() << std::endl;
}

void
CSNBListCompressionFix::runEvent() const
{
  double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(dt))
    M_throw() << "A NAN system event time has been found";
#endif
  
  Sim->dSysTime += locdt;
  
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->dynamics.stream(locdt);

  Sim->freestreamAcc += locdt;

  CGNeighbourList& nblist(dynamic_cast<CGNeighbourList&>
			  (*Sim->dynamics.getGlobals()[cellID]));
  
  dout << "Rebuilding the neighbour list named " << nblist.getName()
	   << "\nNColl = " << Sim->eventCount
	   << "\nSys t = " << Sim->dSysTime / Sim->dynamics.units().unitTime() << std::endl;
  
  nblist.reinitialise(1.0001 * nblist.getMaxSupportedInteractionLength());
  
  dt = (nblist.getMaxSupportedInteractionLength()
	/ nblist.getMaxInteractionLength() - 1.0) / growthRate - Sim->dSysTime;
}


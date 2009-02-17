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

#include "nblistCompressionFix.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/neighbourList.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../../schedulers/scheduler.hpp"

CSNBListCompressionFix::CSNBListCompressionFix(DYNAMO::SimData* nSim, Iflt nGR, size_t nblistID):
  CSystem(nSim),
  growthRate(nGR),
  cellID(nblistID)
{
  sysName = "GlobalCellsCompressionHack";
  type = NON_EVENT;

  if (dynamic_cast<const CGNeighbourList*>(Sim->Dynamics.getGlobals()[cellID].get_ptr()) == NULL)
    D_throw() << "The ID passed to CSNBListCompressionFix isn't a CGNeighbourList";
}

void
CSNBListCompressionFix::initialise(size_t nID)
{
  ID = nID;

  if (dynamic_cast<const CGNeighbourList*>(Sim->Dynamics.getGlobals()[cellID].get_ptr()) == NULL)
    D_throw() << "Have the globals been shuffled? The cellID is no longer a CGNeighbourList.";
  
  const CGNeighbourList& nblist(dynamic_cast<const CGNeighbourList&>
			       (*Sim->Dynamics.getGlobals()[cellID]));

  dt = (nblist.getMaxSupportedInteractionLength() / nblist.getMaxInteractionLength() - 1.0) / growthRate;

  I_cout() << "Compression Hack Loaded"
	   << "\nFor global " << nblist.getName()
	   << "\nCompression rate = " 
	   << growthRate / Sim->Dynamics.units().unitTime()
	   << "\nSim Units compression rate = " << growthRate
	   << "\nMax length of interaction = " 
	   << nblist.getMaxSupportedInteractionLength() / Sim->Dynamics.units().unitLength()
	   << "\nMaximum supported length = "
	   << nblist.getMaxSupportedInteractionLength() / Sim->Dynamics.units().unitLength()
	   << "\nFirst halt scheduled for " 
	   << dt / Sim->Dynamics.units().unitTime();
}

void
CSNBListCompressionFix::runEvent() const
{
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(dt))
    D_throw() << "A NAN system event time has been found";
#endif
  
  Sim->dSysTime += locdt;
  
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->Dynamics.stream(locdt);

  CGNeighbourList& nblist(dynamic_cast<CGNeighbourList&>
			  (*Sim->Dynamics.getGlobals()[cellID]));
  
  I_cout() << "Rebuilding the neighbour list " << nblist.getName()
	   << "\nNColl = " << Sim->lNColl; 
  
  nblist.reinitialise(1.0001 * nblist.getMaxSupportedInteractionLength());
  
  dt = (nblist.getMaxSupportedInteractionLength()
	/ nblist.getMaxInteractionLength() - 1.0) / growthRate;

  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, CNParticleData(), locdt);
}


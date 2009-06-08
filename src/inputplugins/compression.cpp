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

#include "compression.hpp"
#include <boost/foreach.hpp>
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/interactions/squarebond.hpp"
#include "../dynamics/interactions/squarewell.hpp"
#include "../dynamics/interactions/hardsphere.hpp"
#include "../dynamics/ranges/include.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/liouvillean/CompressionL.hpp"
#include "../dynamics/units/units.hpp"
#include "../datatypes/vector.hpp"
#include "../base/is_simdata.hpp"
#include "../schedulers/neighbourlist.hpp"
#include "../dynamics/systems/nblistCompressionFix.hpp"
#include "../dynamics/systems/tHalt.hpp"
#include "../dynamics/species/species.hpp"
#include "../dynamics/globals/gcells.hpp"

CIPCompression::CIPCompression(DYNAMO::SimData* tmp, Iflt GR): 
  CInputPlugin(tmp, "CompressionPlugin"),
  growthRate(GR)
{
  I_cout() << "Compression plugin loaded\n"
	   << "Compaction parameter gamma " << growthRate;
}

void 
CIPCompression::MakeGrowth()
{
  I_cout() << "Backing up old liouvillean";

  //Required to reset the dynamics
  Sim->Dynamics.Liouvillean().updateAllParticles();

  oldLio = Sim->Dynamics.Liouvillean().Clone();

  I_cout() << "Loading compression liouvillean";
  Sim->Dynamics.setLiouvillean(new CLCompression(Sim, growthRate 
						 / (Sim->Dynamics.units()
						    .unitTime())));
}

void
CIPCompression::RestoreSystem()
{
  I_cout() << "Restoring original liouvillean";

  //Required to finish off the compression dynamics
  Sim->Dynamics.Liouvillean().updateAllParticles();

  if (dynamic_cast<CSNeighbourList*>(Sim->ptrScheduler) != NULL)
    {
      BOOST_FOREACH(smrtPlugPtr<CGlobal>& ptr, Sim->Dynamics.getGlobals())
	if (dynamic_cast<const CGCells*>(ptr.get_ptr()) != NULL)      
	  {
	    //Rebulid the collision scheduler without the overlapping cells!
	    CGCells& cells(dynamic_cast<CGCells&>(*ptr));
	    
	    cells.setLambda(oldLambda);
	  }
    }
  else
    I_cout() << "No cellular device to fix";

  Sim->Dynamics.rescaleLengths(Sim->dSysTime * growthRate
			       / Sim->Dynamics.units().unitTime());

  Sim->Dynamics.setLiouvillean(oldLio->Clone());

  Iflt volume = 0.0;
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp->getIntPtr()->hardCoreDiam(), NDIM) * sp->getCount();
  
  Sim->ssHistory << "\nCompression Dynamics run"
		 << "\nEnd packing fraction" 
		 << PI * volume / (6 * Sim->Dynamics.units().simVolume());
}

void
CIPCompression::checkOverlaps()
{
  //Just check the bonds are still valid
  Sim->Dynamics.SystemOverlapTest();
}


void
CIPCompression::CellSchedulerHack()
{
  for (size_t i(0); i < Sim->Dynamics.getGlobals().size(); ++i)
    {      
      if (dynamic_cast<const CGNeighbourList*>(Sim->Dynamics.getGlobals()[i].get_ptr()) != NULL)
	{
	  
	  if (dynamic_cast<const CGCells*>(Sim->Dynamics.getGlobals()[i]
					   .get_ptr()) != NULL)
	    {
	      //Rebulid the collision scheduler without the overlapping
	      //cells as it reduces the number of reinitialisations
	      //required
	      CGCells& cells(static_cast<CGCells&>
			     (*Sim->Dynamics.getGlobals()[i]));	    
	      oldLambda = cells.getLambda();	    
	      cells.setLambda(0.0);
	    }
	  
	  //Add the system watcher
	  Sim->Dynamics.addSystem
	    (new CSNBListCompressionFix(Sim, growthRate 
					/ Sim->Dynamics.units().unitTime(),
					i));
	}
    }
}

void 
CIPCompression::limitPackingFraction(Iflt targetp)
{
  I_cout() << "Limiting maximum packing fraction to " << targetp;
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp->getIntPtr()->hardCoreDiam(), NDIM) * sp->getCount();
  
  Iflt packfrac = PI * volume / (6 * (Sim->Dynamics.units().simVolume()));
  
  if (targetp < packfrac)
    D_throw() << "Target packing fraction is lower than current!";
  
  Sim->Dynamics.addSystem(new CStHalt(Sim, (pow(targetp / packfrac, 1.0/3.0) 
					    - 1.0) / growthRate, 
				      "CompresionLimiter"));
}

void 
CIPCompression::limitDensity(Iflt targetrho)
{
  I_cout() << "Limiting maximum density to " << targetrho;

  //Get the avg molecular volume
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->Dynamics.getSpecies())
    volume += std::pow(sp->getIntPtr()->hardCoreDiam(), static_cast<int>(NDIM)) * sp->getCount();
  
  Iflt molVol = PI * volume / (6.0 * Sim->vParticleList.size()
			       * Sim->Dynamics.units().unitVolume());

  I_cout() << "Corresponding packing fraction for that density is "
	   << molVol * targetrho;
  limitPackingFraction(molVol * targetrho);
}

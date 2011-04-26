/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include "../dynamics/globals/neighbourList.hpp"

CIPCompression::CIPCompression(DYNAMO::SimData* tmp, double GR): 
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
  Sim->dynamics.getLiouvillean().updateAllParticles();

  oldLio = Sim->dynamics.getLiouvillean().Clone();

  I_cout() << "Loading compression liouvillean";
  Sim->dynamics.setLiouvillean(new LCompression(Sim, growthRate 
						 / (Sim->dynamics.units()
						    .unitTime())));
}

void
CIPCompression::RestoreSystem()
{
  I_cout() << "Restoring original liouvillean";

  //Required to finish off the compression dynamics
  Sim->dynamics.getLiouvillean().updateAllParticles();

  if (dynamic_cast<CSNeighbourList*>(Sim->ptrScheduler) != NULL)
    {
      BOOST_FOREACH(magnet::ClonePtr<Global>& ptr, Sim->dynamics.getGlobals())
	if (dynamic_cast<const CGNeighbourList*>(ptr.get_ptr()) != NULL)      
	  //Rebulid the collision scheduler without the overlapping cells!
	  dynamic_cast<CGNeighbourList&>(*ptr).setCellOverlap(true);
    }
  else
    I_cout() << "No cellular device to fix";

  double rescale_factor = 1.0 + Sim->dSysTime * growthRate / Sim->dynamics.units().unitTime();
  Sim->dynamics.units().rescaleLength(rescale_factor);
  Sim->_properties.rescaleUnit(Property::Units::L, rescale_factor);

  Sim->dynamics.setLiouvillean(oldLio->Clone());
  
  Sim->ssHistory << "\nCompression dynamics run"
		 << "\nEnd packing fraction" 
		 << Sim->dynamics.getPackingFraction();
}

void
CIPCompression::checkOverlaps()
{
  //Just check the bonds are still valid
  Sim->dynamics.SystemOverlapTest();
}


void
CIPCompression::CellSchedulerHack()
{
  for (size_t i(0); i < Sim->dynamics.getGlobals().size(); ++i)
    {      
      if (dynamic_cast<const CGNeighbourList*>(Sim->dynamics.getGlobals()[i].get_ptr()) != NULL)
	{
	  //Rebulid the collision scheduler without the overlapping
	  //cells, otherwise cells are always rebuilt as they overlap
	  //such that the maximum supported interaction distance is
	  //equal to the current maximum interaction distance.
	  static_cast<CGNeighbourList&>(*Sim->dynamics.getGlobals()[i]).setCellOverlap(false);
	  
	  //Add the system watcher
	  Sim->dynamics.addSystem
	    (new CSNBListCompressionFix(Sim, growthRate 
					/ Sim->dynamics.units().unitTime(),
					i));
	}
    }
}

void 
CIPCompression::limitPackingFraction(double targetp)
{
  I_cout() << "Limiting maximum packing fraction to " << targetp;
  
  double packfrac = Sim->dynamics.getPackingFraction();
  
  if (targetp < packfrac)
    M_throw() << "Target packing fraction is lower than current!";
  
  Sim->dynamics.addSystem(new CStHalt(Sim, (pow(targetp / packfrac, 1.0/3.0) 
					    - 1.0) / growthRate, 
				      "CompresionLimiter"));
}

void 
CIPCompression::limitDensity(double targetrho)
{
  I_cout() << "Limiting maximum density to " << targetrho;
  
  double molVol = (Sim->dynamics.getPackingFraction() 
		   * Sim->dynamics.units().simVolume())
    / (Sim->N * Sim->dynamics.units().unitVolume());

  I_cout() << "Corresponding packing fraction for that density is "
	   << molVol * targetrho;
  limitPackingFraction(molVol * targetrho);
}

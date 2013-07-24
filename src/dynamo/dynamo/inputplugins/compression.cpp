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

#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/particle.hpp>

#include <dynamo/interactions/squarebond.hpp>
#include <dynamo/interactions/squarewell.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/neighbourlist.hpp>
#include <dynamo/systems/nblistCompressionFix.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/globals/neighbourList.hpp>

namespace dynamo {

  IPCompression::IPCompression(dynamo::Simulation* tmp, double GR): 
    InputPlugin(tmp, "CompressionPlugin"),
    growthRate(GR)
  {
    dout << "Compression plugin loaded\n"
	 << "Compaction parameter gamma " << growthRate << std::endl;
  }

  void 
  IPCompression::MakeGrowth()
  {
    dout << "Backing up old dynamics" << std::endl;

    //Required to reset the dynamics
    Sim->dynamics->updateAllParticles();

    oldLio = Sim->dynamics;

    dout << "Loading compression dynamics" << std::endl;
    Sim->dynamics 
      = shared_ptr<dynamo::Dynamics>
      (new DynCompression(Sim, growthRate / Sim->units.unitTime()));
  }

  void
  IPCompression::RestoreSystem()
  {
    dout << "Restoring original dynamics" << std::endl;

    //Required to finish off the compression dynamics
    Sim->dynamics->updateAllParticles();

    if (std::dynamic_pointer_cast<SNeighbourList>(Sim->ptrScheduler))
      {
	for (shared_ptr<System>& ptr : Sim->systems)
	  if (std::dynamic_pointer_cast<SysNBListCompressionFix>(ptr))
	    static_cast<SysNBListCompressionFix&>(*ptr).fixNBlistForOutput();

	for (shared_ptr<Global>& ptr : Sim->globals)
	  if (std::dynamic_pointer_cast<GNeighbourList>(ptr))
	    //Rebulid the collision scheduler without the overlapping cells!
	    static_cast<GNeighbourList&>(*ptr).setCellOverlap(true);
      }
    else
      dout << "No cellular device to fix" << std::endl;

    double rescale_factor = 1.0 + Sim->systemTime * growthRate / Sim->units.unitTime();

    // The length scale is rescaled as the particles have grown. We want
    // that if a particle had a radius of 1 before the compression, it
    // will have a radius of 1 after the compression (but the simulation
    // volume will be less).
    Sim->units.rescaleLength(rescale_factor);
    // The time scale is also rescaled, so that the energy and velocity
    // scales are unchanged.
    Sim->units.rescaleTime(rescale_factor);
    Sim->_properties.rescaleUnit(Property::Units::L, rescale_factor);
    Sim->_properties.rescaleUnit(Property::Units::T, rescale_factor);

    Sim->dynamics = oldLio;
  }

  void
  IPCompression::CellSchedulerHack()
  {
    for (size_t i(0); i < Sim->globals.size(); ++i)
      {      
	if (std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[i]))
	  {
	    //Rebulid the collision scheduler without the overlapping
	    //cells, otherwise cells are always rebuilt as they overlap
	    //such that the maximum supported interaction distance is
	    //equal to the current maximum interaction distance.
	    static_cast<GNeighbourList&>(*Sim->globals[i]).setCellOverlap(false);
	    
	    //Add the system watcher
	    Sim->systems.push_back
	      (shared_ptr<System>
	       (new SysNBListCompressionFix(Sim, growthRate 
					   / Sim->units.unitTime(),
					   i)));
	  }
      }
  }

  void 
  IPCompression::limitPackingFraction(double targetp)
  {
    dout << "Limiting maximum packing fraction to " << targetp << std::endl;
  
    double packfrac = Sim->getPackingFraction();
  
    if (targetp < packfrac)
      M_throw() << "Target packing fraction is lower than current!";
  
    Sim->systems.push_back
      (shared_ptr<System>
       (new SystHalt(Sim, (pow(targetp / packfrac, 1.0/3.0) - 1.0) / growthRate, 
		    "CompresionLimiter")));
  }

  void 
  IPCompression::limitDensity(double targetrho)
  {
    dout << "Limiting maximum density to " << targetrho << std::endl;
  
    double molVol = (Sim->getPackingFraction() * Sim->getSimVolume())
      / (Sim->N * Sim->units.unitVolume());

    dout << "Corresponding packing fraction for that density is "
	 << molVol * targetrho << std::endl;
    limitPackingFraction(molVol * targetrho);
  }
}

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

#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/interactions/squarebond.hpp>
#include <dynamo/dynamics/interactions/squarewell.hpp>
#include <dynamo/dynamics/interactions/hardsphere.hpp>
#include <dynamo/dynamics/ranges/include.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/liouvillean/CompressionL.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/schedulers/neighbourlist.hpp>
#include <dynamo/dynamics/systems/nblistCompressionFix.hpp>
#include <dynamo/dynamics/systems/tHalt.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <boost/foreach.hpp>

namespace dynamo {

  CIPCompression::CIPCompression(dynamo::SimData* tmp, double GR): 
    CInputPlugin(tmp, "CompressionPlugin"),
    growthRate(GR)
  {
    dout << "Compression plugin loaded\n"
	 << "Compaction parameter gamma " << growthRate << std::endl;
  }

  void 
  CIPCompression::MakeGrowth()
  {
    dout << "Backing up old liouvillean" << std::endl;

    //Required to reset the dynamics
    Sim->dynamics.getLiouvillean().updateAllParticles();

    oldLio = Sim->dynamics.getLiouvilleanPtr();

    dout << "Loading compression liouvillean" << std::endl;
    Sim->dynamics.setLiouvillean(std::tr1::shared_ptr<dynamo::Liouvillean>
				 (new LCompression(Sim, growthRate 
						   / Sim->dynamics.units().unitTime())));
  }

  void
  CIPCompression::RestoreSystem()
  {
    dout << "Restoring original liouvillean" << std::endl;

    //Required to finish off the compression dynamics
    Sim->dynamics.getLiouvillean().updateAllParticles();

    if (std::tr1::dynamic_pointer_cast<SNeighbourList>(Sim->ptrScheduler))
      {
	BOOST_FOREACH(std::tr1::shared_ptr<Global>& ptr, Sim->dynamics.getGlobals())
	  if (std::tr1::dynamic_pointer_cast<GNeighbourList>(ptr))
	    //Rebulid the collision scheduler without the overlapping cells!
	    static_cast<GNeighbourList&>(*ptr).setCellOverlap(true);
      }
    else
      dout << "No cellular device to fix" << std::endl;

    double rescale_factor = 1.0 + Sim->dSysTime * growthRate / Sim->dynamics.units().unitTime();

    // The length scale is rescaled as the particles have grown. We want
    // that if a particle had a radius of 1 before the compression, it
    // will have a radius of 1 after the compression (but the simulation
    // volume will be less).
    Sim->dynamics.units().rescaleLength(rescale_factor);
    // The time scale is also rescaled, so that the energy and velocity
    // scales are unchanged.
    Sim->dynamics.units().rescaleTime(rescale_factor);
    Sim->_properties.rescaleUnit(Property::Units::L, rescale_factor);
    Sim->_properties.rescaleUnit(Property::Units::T, rescale_factor);

    Sim->dynamics.getLiouvilleanPtr() = oldLio;
  
    Sim->ssHistory << "\nCompression dynamics run"
		   << "\nEnd packing fraction" 
		   << Sim->dynamics.getPackingFraction();
  }

  void
  CIPCompression::CellSchedulerHack()
  {
    for (size_t i(0); i < Sim->dynamics.getGlobals().size(); ++i)
      {      
	if (std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[i]))
	  {
	    //Rebulid the collision scheduler without the overlapping
	    //cells, otherwise cells are always rebuilt as they overlap
	    //such that the maximum supported interaction distance is
	    //equal to the current maximum interaction distance.
	    static_cast<GNeighbourList&>(*Sim->dynamics.getGlobals()[i]).setCellOverlap(false);
	    
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
    dout << "Limiting maximum packing fraction to " << targetp << std::endl;
  
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
    dout << "Limiting maximum density to " << targetrho << std::endl;
  
    double molVol = (Sim->dynamics.getPackingFraction() * Sim->dynamics.getSimVolume())
      / (Sim->N * Sim->dynamics.units().unitVolume());

    dout << "Corresponding packing fraction for that density is "
	 << molVol * targetrho << std::endl;
    limitPackingFraction(molVol * targetrho);
  }
}

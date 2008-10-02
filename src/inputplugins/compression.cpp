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
#include "../schedulers/cells.hpp"
#include "../dynamics/systems/compressionhack.hpp"
#include "../dynamics/systems/tHalt.hpp"
#include "../dynamics/species/species.hpp"

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

  {
    CSCells* tmpPtr = dynamic_cast<CSCells*>(Sim->ptrScheduler);
    if (tmpPtr != NULL)
      tmpPtr->setLambda(oldLambda);
  }

  Sim->Dynamics.rescaleLengths(Sim->dSysTime * growthRate 
			       / Sim->Dynamics.units().unitTime());

  Sim->Dynamics.setLiouvillean(oldLio->Clone());

  Iflt volume = 0.0;
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp.getIntPtr()->hardCoreDiam(), NDIM) * sp.getCount();
  
  Sim->ssHistory << "\nCompression Dynamics run"
		 << "\nEnd packing fraction" << PI * volume / (6 * Sim->Dynamics.units().simVolume());
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
  CSCells* tmpPtr = dynamic_cast<CSCells*>(Sim->ptrScheduler);
  if (tmpPtr == NULL)
    {
      I_cout() << "Cannot implement cellular compression hack, not cellular";
      return;
    }
  
  //Rebulid the collision scheduler without the overlapping cells!
  oldLambda = tmpPtr->getLambda();
  tmpPtr->setLambda(0.0);
  //Add the system watcher
  Sim->Dynamics.addSystem(new CSCellHack(Sim, growthRate 
					 / Sim->Dynamics.units().unitTime()));
}

void 
CIPCompression::limitPackingFraction(Iflt targetp)
{
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp.getIntPtr()->hardCoreDiam(), NDIM) * sp.getCount();
  
  Iflt packfrac = PI * volume / (6 * (Sim->Dynamics.units().simVolume()));

  if (targetp < packfrac)
    I_throw() << "Target packing fraction is lower than current!";
  
  Sim->Dynamics.addSystem(new CStHalt(Sim, (pow(targetp / packfrac, 1.0/3.0) 
					    - 1.0) / growthRate, 
				      "CompresionLimiter"));
}

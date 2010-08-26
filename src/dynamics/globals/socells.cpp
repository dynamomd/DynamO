/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "socells.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include <boost/static_assert.hpp>
#include <boost/math/special_functions/pow.hpp>
#include "../liouvillean/NewtonianGravityL.hpp"


CGSOCells::CGSOCells(DYNAMO::SimData* nSim, const std::string& name):
  Global(nSim, "SingleOccupancyCells"),
  cellCount(0),
  cellDimension(1,1,1),
  cuberootN(0)
{
  globName = name;
  I_cout() << "Single occupancy cells loaded";
}

CGSOCells::CGSOCells(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  Global(ptrSim, "SingleOccupancyCells"),
  cellCount(0),
  cellDimension(1,1,1),
  cuberootN(0)
{
  operator<<(XML);

  I_cout() << "Single occupancy cells loaded";
}

void 
CGSOCells::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGSOCells";
    }
}

GlobalEvent 
CGSOCells::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  //This 
  //Sim->dynamics.getLiouvillean().updateParticle(part);
  //is not required as we compensate for the delay using 
  //Sim->dynamics.getLiouvillean().getParticleDelay(part)

  Vector CellOrigin;
  size_t ID(part.getID());

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->aspectRatio[iDim];
      ID /= cuberootN;
    }

  return GlobalEvent(part,
		    Sim->dynamics.getLiouvillean().
		    getSquareCellCollision2
		    (part, CellOrigin,
		     cellDimension)
		    -Sim->dynamics.getLiouvillean().getParticleDelay(part),
		    CELL, *this);
}

void
CGSOCells::runEvent(const Particle& part) const
{
  Sim->dynamics.getLiouvillean().updateParticle(part);

  Vector CellOrigin;
  size_t ID(part.getID());

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->aspectRatio[iDim];
      ID /= cuberootN;
    }
  
  //Determine the cell transition direction, its saved
  int cellDirectionInt(Sim->dynamics.getLiouvillean().
		       getSquareCellCollision3
		       (part, CellOrigin, 
			cellDimension));

  size_t cellDirection = abs(cellDirectionInt) - 1;

  GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif

  Sim->dSysTime += iEvent.getdt();
    
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  Sim->dynamics.stream(iEvent.getdt());

  Vector vNorm(0,0,0);

  Vector pos(part.getPosition()), vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(pos, vel);

  vNorm[cellDirection] = (cellDirectionInt > 0) ? -1 : +1; 
    
  //Run the collision and catch the data
  NEventData EDat(Sim->dynamics.getLiouvillean().runWallCollision
		      (part, vNorm, 1.0));

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

}

void 
CGSOCells::initialise(size_t nID)
{
  ID=nID;
  
  cuberootN = (unsigned long)(std::pow(Sim->lN, 1.0/3.0) + 0.5);
  
  if (boost::math::pow<3>(cuberootN) != Sim->lN)
    D_throw() << "Cannot use single occupancy cells without a integer cube root of N"
	      << "\nN = " << Sim->lN
	      << "\nN^(1/3) = " << cuberootN;

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    cellDimension[iDim] = Sim->aspectRatio[iDim] / cuberootN;

  if (Sim->dynamics.liouvilleanTypeTest<LNewtonianGravity>())
    I_cout() << "Warning, in order for SingleOccupancyCells to work in gravity\n"
	     << "You must add the ParabolaSentinel Global event.";

}

void
CGSOCells::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SOCells"
      << xmlw::attr("Name") << globName;
}

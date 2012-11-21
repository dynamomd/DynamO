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

#include <dynamo/globals/socells.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/particle.hpp>
#include <boost/static_assert.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  GSOCells::GSOCells(dynamo::Simulation* nSim, const std::string& name):
    Global(nSim, "SingleOccupancyCells"),
    cellDimension(1,1,1),
    cuberootN(0)
  {
    globName = name;
    dout << "Single occupancy cells loaded" << std::endl;
  }

  GSOCells::GSOCells(const magnet::xml::Node&XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "SingleOccupancyCells"),
    cellDimension(1,1,1),
    cuberootN(0)
  {
    operator<<(XML);

    dout << "Single occupancy cells loaded" << std::endl;
  }

  void 
  GSOCells::operator<<(const magnet::xml::Node& XML)
  {
    globName = XML.getAttribute("Name");	
  }

  GlobalEvent 
  GSOCells::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    //This 
    //Sim->dynamics->updateParticle(part);
    //is not required as we compensate for the delay using 
    //Sim->dynamics->getParticleDelay(part)

    Vector CellOrigin;
    size_t ID(part.getID());

    for (size_t iDim(0); iDim < NDIM; ++iDim)
      {
	CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->primaryCellSize[iDim];
	ID /= cuberootN;
      }

    return GlobalEvent(part,
		       Sim->dynamics->
		       getSquareCellCollision2
		       (part, CellOrigin,
			cellDimension)
		       -Sim->dynamics->getParticleDelay(part),
		       CELL, *this);
  }

  void
  GSOCells::runEvent(Particle& part, const double) const
  {
    Sim->dynamics->updateParticle(part);

    Vector CellOrigin;
    size_t ID(part.getID());

    for (size_t iDim(0); iDim < NDIM; ++iDim)
      {
	CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->primaryCellSize[iDim];
	ID /= cuberootN;
      }
  
    //Determine the cell transition direction, its saved
    int cellDirectionInt(Sim->dynamics->
			 getSquareCellCollision3
			 (part, CellOrigin, 
			  cellDimension));

    size_t cellDirection = abs(cellDirectionInt) - 1;

    GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(iEvent.getdt()))
      M_throw() << "A NAN Interaction collision time has been found"
		<< iEvent.stringData(Sim);
  
    if (iEvent.getdt() == HUGE_VAL)
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		<< iEvent.stringData(Sim);
#endif

    Sim->systemTime += iEvent.getdt();
    
    Sim->ptrScheduler->stream(iEvent.getdt());
  
    Sim->stream(iEvent.getdt());

    Vector vNorm(0,0,0);

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    vNorm[cellDirection] = (cellDirectionInt > 0) ? -1 : +1; 
    
    //Run the collision and catch the data
    NEventData EDat(Sim->dynamics->runPlaneEvent(part, vNorm, 1.0, 0.0));

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);

  }

  void 
  GSOCells::initialise(size_t nID)
  {
    ID=nID;
  
    cuberootN = (unsigned long)(std::pow(Sim->N, 1.0/3.0) + 0.5);
  
    if (boost::math::pow<3>(cuberootN) != Sim->N)
      M_throw() << "Cannot use single occupancy cells without a integer cube root of N"
		<< "\nN = " << Sim->N
		<< "\nN^(1/3) = " << cuberootN;

    for (size_t iDim(0); iDim < NDIM; ++iDim)
      cellDimension[iDim] = Sim->primaryCellSize[iDim] / cuberootN;

    if (std::tr1::dynamic_pointer_cast<const DynGravity>(Sim->dynamics))
      dout << "Warning, in order for SingleOccupancyCells to work in gravity\n"
	   << "You must add the ParabolaSentinel Global event." << std::endl;

  }

  void
  GSOCells::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "SOCells"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::endtag("Global");
  }
}

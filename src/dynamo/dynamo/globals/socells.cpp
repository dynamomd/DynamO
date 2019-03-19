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
#include <boost/math/special_functions/pow.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  GSOCells::GSOCells(dynamo::Simulation* nSim, const std::string& name):
    Global(nSim, "SingleOccupancyCells"),
    cellDimension({1,1,1})
  {
    globName = name;
    load_cell_origins(nSim->particles);
    dout << "Single occupancy cells loaded" << std::endl;
  }

  GSOCells::GSOCells(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "SingleOccupancyCells"),
    cellDimension({1,1,1})
  {
    operator<<(XML);

    if (cell_origins.empty()) {
      derr << "Loading SOCells from the current particle positions!" << std::endl;
      load_cell_origins(ptrSim->particles);
    }
    
    dout << "Single occupancy cells loaded" << std::endl;
  }

  void 
  GSOCells::operator<<(const magnet::xml::Node& XML)
  {
    globName = XML.getAttribute("Name");
    cellDimension << XML.getNode("CellSize");
    cellDimension *= Sim->units.unitLength();

    if (XML.hasNode("CellOrigins")) {
      Vector pos;
      for (magnet::xml::Node node = XML.getNode("CellOrigins").findNode("Origin"); node.valid(); ++node) {
	pos << node;
	pos *= Sim->units.unitLength();
	cell_origins.push_back(pos);
      }

      if (cell_origins.size() != Sim->N())
	M_throw() << "Number of CellOrigins (" << cell_origins.size() << ") does not match number of particles (" << Sim->N() << ")\n" << XML.getPath();
    }

  }

  Event 
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

    const Vector CellOrigin = cell_origins[part.getID()];
    
    return Event(part, Sim->dynamics->getSquareCellCollision2(part, CellOrigin, cellDimension) - Sim->dynamics->getParticleDelay(part), GLOBAL, CELL, ID);
  }

  void
  GSOCells::runEvent(Particle& part, const double)
  {
    Sim->dynamics->updateParticle(part);

    const Vector CellOrigin = cell_origins[part.getID()];
  
    //Determine the cell transition direction, its saved
    int cellDirectionInt(Sim->dynamics->getSquareCellCollision3(part, CellOrigin, cellDimension));

    size_t cellDirection = abs(cellDirectionInt) - 1;

    Event iEvent = getEvent(part);

#ifdef DYNAMO_DEBUG 
    if (std::isnan(iEvent._dt))
      M_throw() << "A NAN Interaction collision time has been found";
  
    if (iEvent._dt == std::numeric_limits<float>::infinity())
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n";
#endif

    Sim->systemTime += iEvent._dt;
    
    Sim->ptrScheduler->stream(iEvent._dt);
  
    Sim->stream(iEvent._dt);

    Vector vNorm({0,0,0});

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    vNorm[cellDirection] = (cellDirectionInt > 0) ? -1 : +1; 
    
    //Run the collision and catch the data
    NEventData EDat(Sim->dynamics->runPlaneEvent(part, vNorm, 1.0, 0.0));

    Sim->_sigParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  void 
  GSOCells::initialise(size_t nID)
  {
    Global::initialise(nID);

    if (std::dynamic_pointer_cast<const DynGravity>(Sim->dynamics))
      dout << "Warning, in order for SingleOccupancyCells to work in gravity\n"
	   << "You must add the ParabolaSentinel Global event." << std::endl;
  }

  void
  GSOCells::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "SOCells"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::tag("CellSize") << cellDimension / Sim->units.unitLength()
	<< magnet::xml::endtag("CellSize");

    XML << magnet::xml::tag("CellOrigins");
    for (const Vector cellorigin : cell_origins)
      XML << magnet::xml::tag("Origin")
	  << cellorigin / Sim->units.unitLength()
	  << magnet::xml::endtag("Origin");
    XML << magnet::xml::endtag("CellOrigins");
      
    
    XML << magnet::xml::endtag("Global");
  }

  
  void
  GSOCells::load_cell_origins(const std::vector<Particle> particles) {
    cell_origins.resize(Sim->particles.size());
    
    for (const Particle& p : particles) {
      cell_origins[p.getID()] = p.getPosition() - 0.5 * cellDimension;
    }
  }
}




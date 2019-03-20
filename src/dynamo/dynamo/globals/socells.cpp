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
    _cellD(0)
  {
    globName = name;
    load_cell_origins(nSim->particles);
    dout << "Single occupancy cells loaded" << std::endl;
  }

  GSOCells::GSOCells(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "SingleOccupancyCells"),
    _cellD(0)
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

    if (XML.hasAttribute("Diameter"))
      _cellD = XML.getAttribute("Diameter").as<double>() * Sim->units.unitLength();
  }

  Event 
  GSOCells::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    //Not needed as we compensate for particle delay
    //Sim->dynamics->updateParticle(part);

    //Create a fake particle which represents the cell center
    const Particle cellParticle(cell_origins[part.getID()], Vector({0,0,0}), -1);

    Vector pos = part.getPosition() - cellParticle.getPosition();
    Sim->BCs->applyBC(pos); //We don't apply the PBC, as 

    //if (part.getID() == 2) {
    //  dout << "#Testing event, Particle " << part.getID() << std::endl;
    //  dout << "#Position " << pos / Sim->units.unitLength() << std::endl;
    //  dout << "#Velocity " << part.getVelocity() / Sim->units.unitVelocity() << std::endl;
    //  dout << "#Distance " << pos.nrm() / Sim->units.unitLength() << std::endl;
    //  dout << "#Delay " << Sim->dynamics->getParticleDelay(part) / Sim->units.unitTime() << std::endl;
    //  double dt = Sim->dynamics->SphereSphereOutRoot(part, cellParticle, _cellD);
    //  dout << "#event t = " << dt / Sim->units.unitTime() << std::endl;
    //  Vector finalP = (pos + dt * part.getVelocity());
    //  dout << "#Final Position " << finalP / Sim->units.unitLength() << std::endl;
    //  dout << "#Relative final p" << (finalP.nrm() - _cellD)/_cellD<< std::endl;
    //  dout << " cell origin " << cell_origins[part.getID()] / Sim->units.unitLength() << std::endl;
    //}
        
    return Event(part, Sim->dynamics->SphereSphereOutRoot(part, cellParticle, _cellD) - Sim->dynamics->getParticleDelay(part), GLOBAL, CELL, ID);
  }

  void
  GSOCells::runEvent(Particle& part, const double)
  {    
    Sim->dynamics->updateParticle(part);
    Sim->ptrScheduler->popNextEvent();
    Event iEvent = getEvent(part);

#ifdef DYNAMO_DEBUG 
    if (std::isnan(iEvent._dt))
      M_throw() << "A NAN Interaction collision time has been found";
  
    if (iEvent._dt == std::numeric_limits<float>::infinity())
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n";
#endif
    
    //Move the system forward to the time of the event
    Sim->systemTime += iEvent._dt;
    Sim->ptrScheduler->stream(iEvent._dt);
    Sim->stream(iEvent._dt);
    ++Sim->eventCount;

    Vector pos = part.getPosition() - cell_origins[part.getID()];
    Sim->BCs->applyBC(pos); //We don't apply the PBC, as 
    
    //dout << "!Particle " << part.getID() << " at " << part.getPosition() / Sim->units.unitLength() << std::endl;
    //dout << "!Cell origin " << cell_origins[part.getID()] / Sim->units.unitLength() << std::endl;
    //dout << "!Normal " << pos.normal() / Sim->units.unitLength() << std::endl;
    //dout << "!Distance " << pos.nrm() / Sim->units.unitLength() << std::endl;
    //dout << "!CellD " << _cellD / Sim->units.unitLength() << std::endl;
    //dout << "!Relative distance " << pos.nrm() / _cellD << std::endl;
    //dout << "!Perp velocity " << (pos.normal() | part.getVelocity()) << std::endl;

    //Now execute the event
    NEventData EDat(Sim->dynamics->runPlaneEvent(part, pos.normal(), 1.0, _cellD));
    //dout << "!Perp velocity post " << (pos.normal() | part.getVelocity()) << std::endl;

    if (pos.nrm() > _cellD * 1.00000001)
      derr << "Particle " << part.getID() << " outside the cell by " << (pos.nrm() - _cellD) / _cellD  << std::endl;
    if (pos.nrm() < _cellD * 0.99999999)
      derr << "Particle " << part.getID() << " inside the cell by " << (pos.nrm() - _cellD) / _cellD  << std::endl;
    
    //Now we're past the event update everything
    Sim->_sigParticleUpdate(EDat);
    Sim->ptrScheduler->fullUpdate(part);
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);

  }

  void 
  GSOCells::initialise(size_t nID)
  {
    Global::initialise(nID);

    //If not set already, set the diameter of the cells such that the
    //simulation volume and total cell volume are equal.
    if (_cellD == 0) {
      const double cellVolume = Sim->getSimVolume() / Sim->N();
      _cellD = std::cbrt(cellVolume * 6 / M_PI);
    }

    if (((_cellD >= 0.5 * Sim->primaryCellSize[0])
	|| (_cellD >= 0.5 * Sim->primaryCellSize[1])
	|| (_cellD >= 0.5 * Sim->primaryCellSize[2]))
	&& (std::dynamic_pointer_cast<BCPeriodic>(Sim->BCs)))
      M_throw() << "ERROR: SOCells diameter (" << _cellD / Sim->units.unitLength() << ") is more than half the primary image size (" << Sim->primaryCellSize << "), this will break in periodic boundary conditions";

    for (const Particle& p : Sim->particles) {
      Vector pos = p.getPosition() - cell_origins[p.getID()];
      Sim->BCs->applyBC(pos);
      if (pos.nrm2() > _cellD * _cellD)
	derr << "Particle " << p.getID() << " is at a distance of "
	     << (pos).nrm() / Sim->units.unitLength()
	     << " cell origin " << cell_origins[p.getID()] / Sim->units.unitLength()
	     << " outside its SOCell where the diameter is " << _cellD / Sim->units.unitLength() << std::endl;
    }
  }

  void
  GSOCells::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "SOCells"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("Diameter") << _cellD / Sim->units.unitLength();

    XML << magnet::xml::tag("CellOrigins");
    for (const Vector cellorigin : cell_origins)
      XML << magnet::xml::tag("Origin")
	  << cellorigin / Sim->units.unitLength()
	  << magnet::xml::endtag("Origin");
    XML << magnet::xml::endtag("CellOrigins")
	<< magnet::xml::endtag("Global");
  }

  
  void
  GSOCells::load_cell_origins(const std::vector<Particle> particles) {
    cell_origins.resize(Sim->particles.size());
    
    for (const Particle& p : particles) {
      cell_origins[p.getID()] = p.getPosition();
    }
  }
}




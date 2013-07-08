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

#ifdef DYNAMO_visualizer

#include <dynamo/interactions/glyphrepresentation.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/systems/visualizer.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/coilRenderObj.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/include.hpp>
#include <coil/clWindow.hpp>
#include <coil/RenderObj/DataSet.hpp>
#include <algorithm>

namespace dynamo {
  SVisualizer::SVisualizer(dynamo::Simulation* nSim, std::string nName, double tickFreq):
    System(nSim)
  {
    //Convert to output units of time
    tickFreq /= Sim->units.unitTime();
    //Stop zero tick times, replace them with a guess of 0.01
    tickFreq += (tickFreq == 0) * 0.01;

    //We want to ensure we get at least one update before anything
    //occurs in the system.
    dt = -HUGE_VAL;
  
    sysName = "Visualizer";

    //Build a window, ready to display it
    _window.reset(new coil::CLGLWindow("Visualizer : " + nName, tickFreq, true));

    //Add all of the objects to be rendered
    _particleData.reset(new coil::DataSet("Particles", Sim->N, GlyphRepresentation::SPHERE_GLYPH));
    _window->addRenderObj(_particleData);

    for (shared_ptr<Local>& local : Sim->locals)
      {
	CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));
	if (obj != NULL) _window->addRenderObj(obj->getCoilRenderObj());
      }

    //Initialise coil
    _coil.getInstance().addWindow(_window);

    //Now initialise the data sets
    initDataSet();
    _window->_updateDataSignal.connect<SVisualizer, &SVisualizer::updateRenderData>(this);
    updateRenderData();

    _lastUpdate = boost::posix_time::microsec_clock::local_time();
    /* Now request that the visualiser rescales to the best dimensions for the current system */
    _window->autoscaleView();

    dout << "Visualizer initialised" << std::endl;
  }

  void
  SVisualizer::runEvent() const
  {
    //Dont rewind in time, the -HUGE_VAL time is only used to ensure
    //the event takes place before any event, including negative time events.
    if (dt == -HUGE_VAL) dt = 0;
  
    //Actually move forward the system time
    Sim->systemTime += dt;
    Sim->ptrScheduler->stream(dt);
    //dynamics must be updated first
    Sim->stream(dt);

    if (_window->dynamoParticleSync())
      Sim->dynamics->updateAllParticles();

    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, NEventData(), dt);
  
    std::shared_ptr<DynGravity> dynamics = std::dynamic_pointer_cast<DynGravity>(Sim->dynamics);
    if (dynamics)
      _window->getCamera().setUp(-dynamics->getGravityVector(), true);

    _window->simupdateTick(Sim->systemTime / Sim->units.unitTime());

    //Now that the update has been performed, set up the next "tick"
    dt = _window->getUpdateInterval() * Sim->units.unitTime();
    _lastUpdate = boost::posix_time::microsec_clock::local_time();
  }

  void 
  SVisualizer::initialise(size_t nID)
  { 
    ID = nID;
    (*Sim->_sigParticleUpdate).connect<SVisualizer, &SVisualizer::particlesUpdated>(this);
  }

  void 
  SVisualizer::particlesUpdated(const NEventData&)
  {
    if ((boost::posix_time::microsec_clock::local_time() - _lastUpdate) 
	> boost::posix_time::milliseconds(500))
      {
	dt = -HUGE_VAL;
	Sim->ptrScheduler->rebuildSystemEvents();
      }
  }

  void
  SVisualizer::initDataSet() const
  {
    _particleData->addAttribute("Position", coil::Attribute::COORDINATE | coil::Attribute::DEFAULT_GLYPH_POSITION, 3);
    _particleData->addAttribute("Velocity", coil::Attribute::INTENSIVE, 3);
    _particleData->addAttribute("Size", coil::Attribute::INTENSIVE | coil::Attribute::DEFAULT_GLYPH_SCALING, 3);
    _particleData->addAttribute("Mass", coil::Attribute::EXTENSIVE, 1);
    _particleData->addAttribute("Event Count", coil::Attribute::EXTENSIVE, 1);
    _particleData->addAttribute("ID", coil::Attribute::INTENSIVE | coil::Attribute::DEFAULT_GLYPH_COLOUR, 1);

    _particleData->setPeriodicVectors(Vector(Sim->primaryCellSize[0], 0, 0),
				    Vector(0, Sim->primaryCellSize[1], 0),
				    Vector(0, 0, Sim->primaryCellSize[2]));

    if (Sim->dynamics->hasOrientationData())
      {
	_particleData->addAttribute("Orientation", coil::Attribute::EXTENSIVE | coil::Attribute::DEFAULT_GLYPH_ORIENTATION, 4);
	_particleData->addAttribute("Angular Velocity", coil::Attribute::EXTENSIVE, 3);
      }

    std::vector<GLfloat>& masses = (*_particleData)["Mass"];
    std::vector<GLfloat>& IDs = (*_particleData)["ID"];
    for (const Particle& p : Sim->particles)
      {
	IDs[p.getID()] = p.getID();
	masses[p.getID()] = Sim->species[p]->getMass(p.getID()) / Sim->units.unitMass();
      }
    (*_particleData)["Mass"].flagNewData();
    (*_particleData)["ID"].flagNewData();

    for (auto& species : Sim->species)
      {
	std::vector<GLuint> ids;
	ids.reserve(species->getRange()->size());
	for (auto ID : *species->getRange()) ids.push_back(ID);
	const GlyphRepresentation& data = dynamic_cast<const GlyphRepresentation&>(*species->getIntPtr());
	_particleData->addPoints(species->getName(), ids, data.getDefaultGlyphType());
      }
  }

  void
  SVisualizer::updateRenderData() const
  {
    if (!_particleData)
      M_throw() << "Updating before the render object has been fetched";
    
    shared_ptr<BCLeesEdwards> BC = std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs);
    if (BC)
      _particleData->setPeriodicVectors(Vector(Sim->primaryCellSize[0], 0, 0),
				      Vector(BC->getBoundaryDisplacement(), Sim->primaryCellSize[1], 0),
				      Vector(0, 0, Sim->primaryCellSize[2]));


    ///////////////////////POSITION DATA UPDATE
    //Check if the system is compressing and adjust the radius scaling factor
    std::vector<GLfloat>& posdata = (*_particleData)["Position"];
    std::vector<GLfloat>& veldata = (*_particleData)["Velocity"];
    std::vector<GLfloat>& sizes = (*_particleData)["Size"];
    std::vector<GLfloat>& eventCounts = (*_particleData)["Event Count"];
    const std::vector<size_t>& simEventCounts = Sim->ptrScheduler->getEventCounts();
    
    for (const Particle& p : Sim->particles)
      {
	Vector vel = p.getVelocity() / Sim->units.unitVelocity();
	Vector pos = p.getPosition() / Sim->units.unitLength();
	Sim->BCs->applyBC(pos, vel);
	  
	for (size_t i(0); i < NDIM; ++i)
	  {
	    posdata[3 * p.getID() + i] = pos[i];
	    veldata[3 * p.getID() + i] = vel[i];
	  }

	eventCounts[p.getID()] = 0;
	if (!simEventCounts.empty()) 
	  eventCounts[p.getID()] = simEventCounts[p.getID()];
      }

    float rfactor = 1.0;
    if (std::dynamic_pointer_cast<DynCompression>(Sim->dynamics))
      rfactor *= (1 + static_cast<const DynCompression&>(*Sim->dynamics).getGrowthRate() * Sim->systemTime);
    rfactor /= Sim->units.unitLength();

    for (auto& species : Sim->species)
      {
	const GlyphRepresentation& data = dynamic_cast<const GlyphRepresentation&>(*species->getIntPtr());
	for (auto ID : *species->getRange())
	  {
	    Vector psize = rfactor * data.getGlyphSize(ID);
	    for (size_t i(0); i < NDIM; ++i)
	      sizes[3 * ID + i] = psize[i];
	  }
      }

    if (Sim->dynamics->hasOrientationData())
      {
	std::vector<GLfloat>& orientationdata = (*_particleData)["Orientation"];
	std::vector<GLfloat>& angularvdata = (*_particleData)["Angular Velocity"];
	const std::vector<Dynamics::rotData>& data = Sim->dynamics->getCompleteRotData();
	for (const Particle& p : Sim->particles)
	  {
	    for (size_t i(0); i < NDIM; ++i)
	      {
		angularvdata[3 * p.getID() + i] = data[p.getID()].angularVelocity[i] * Sim->units.unitTime();
		orientationdata[4 * p.getID() + i] = data[p.getID()].orientation.imaginary()[i];
	      }
	    orientationdata[4 * p.getID() + 3] = data[p.getID()].orientation.real();
	  }
	(*_particleData)["Angular Velocity"].flagNewData();
	(*_particleData)["Orientation"].flagNewData();
      }

    (*_particleData)["Position"].flagNewData();
    (*_particleData)["Velocity"].flagNewData();
    (*_particleData)["Size"].flagNewData();
    (*_particleData)["Event Count"].flagNewData();
  }
}
#endif

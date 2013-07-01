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
  
    for (const shared_ptr<Species>& spec : Sim->species)
      _window->addRenderObj(spec->createDataSet());

    for (shared_ptr<Local>& local : Sim->locals)
      {
	CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));

	if (obj != NULL)
	  _window->addRenderObj(obj->getCoilRenderObj());
      }

    _coil.getInstance().addWindow(_window);

    for (const shared_ptr<Species>& spec : Sim->species)
      {
	spec->initDataSet();
	_window->_updateDataSignal.connect<Species, &Species::updateRenderData>(spec.get());
	spec->updateRenderData();
      }
  
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
      _window->getCamera().setTransformation(magnet::math::Quaternion::fromToVector(Vector(0,-1,0), dynamics->getGravityVector()));

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
}
#endif

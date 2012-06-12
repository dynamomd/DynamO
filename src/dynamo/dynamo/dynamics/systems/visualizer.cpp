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

#include <dynamo/dynamics/systems/visualizer.hpp>
#include <dynamo/dynamics/liouvillean/CompressionL.hpp>
#include <dynamo/dynamics/coilRenderObj.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/include.hpp>
#include <coil/clWindow.hpp>
#include <coil/RenderObj/DataSet.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

namespace dynamo {
  SVisualizer::SVisualizer(dynamo::SimData* nSim, std::string nName, double tickFreq):
    System(nSim)
  {
    //Convert to output units of time
    tickFreq /= Sim->dynamics.units().unitTime();
    //Stop zero tick times 
    tickFreq += (tickFreq == 0);

    //We want to ensure we get at least one update before anything
    //occurs in the system.
    dt = -HUGE_VAL;
  
    sysName = "Visualizer";

    //Build a window, ready to display it
    _window.reset(new coil::CLGLWindow("Visualizer : " + nName, tickFreq, true));
  
    BOOST_FOREACH(const shared_ptr<Species>& spec, Sim->species)
      _window->addRenderObj(spec->createDataSet());

    BOOST_FOREACH(shared_ptr<Local>& local, Sim->dynamics.getLocals())
      {
	CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));

	if (obj != NULL)
	  _window->addRenderObj(obj->getCoilRenderObj());
      }

    _coil.getInstance().addWindow(_window);

    BOOST_FOREACH(const shared_ptr<Species>& spec, Sim->species)
      {
	spec->initDataSet();
	_window->signal_data_update().connect(boost::bind(&Species::updateRenderData, spec.get()));
	spec->updateRenderData();
      }
  
    _lastUpdate = boost::posix_time::microsec_clock::local_time();

    dout << "Visualizer initialised" << std::endl;
  }

  void
  SVisualizer::runEvent() const
  {
    //Dont rewind in time, the -HUGE_VAL time is only used to ensure
    //the event takes place before any event, including negative time events.
    if (dt == -HUGE_VAL) dt = 0;
  
    //Actually move forward the system time
    Sim->dSysTime += dt;
    Sim->ptrScheduler->stream(dt);
    //dynamics must be updated first
    Sim->dynamics.stream(dt);

    double locdt = dt + Sim->freestreamAcc;
    Sim->freestreamAcc = 0;

    if (_window->dynamoParticleSync())
      Sim->liouvillean->updateAllParticles();

    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(*this, NEventData(), locdt);
  
    _window->simupdateTick(Sim->dSysTime / Sim->dynamics.units().unitTime());

    //Now that the update has been performed, set up the next "tick"
    dt = _window->getUpdateInterval() * Sim->dynamics.units().unitTime();
    _lastUpdate = boost::posix_time::microsec_clock::local_time();
  }

  void 
  SVisualizer::initialise(size_t nID)
  { 
    ID = nID; 

    Sim->registerParticleUpdateFunc
      (magnet::function::MakeDelegate(this, &SVisualizer::particlesUpdated));
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

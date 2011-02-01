/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#ifdef DYNAMO_visualizer

#include "visualizer.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../outputplugins/tickerproperty/ticker.hpp"
#include "../units/units.hpp"
#include "../../schedulers/scheduler.hpp"
#include <boost/foreach.hpp>
#include <algorithm>
#include "../../dynamics/include.hpp"
#include <coil/clWindow.hpp>
#include <coil/RenderObj/Function.hpp>
#include <coil/RenderObj/Spheres.hpp>
#include <coil/RenderObj/Lines.hpp>
#include <magnet/CL/CLGL.hpp>
#include "../liouvillean/CompressionL.hpp"

SVisualizer::SVisualizer(DYNAMO::SimData* nSim, std::string nName, double tickFreq):
  System(nSim)
{
  _updateTime = tickFreq * Sim->dynamics.units().unitTime();
  dt = 0;

  sysName = "Visualizer";

  //Build a window, ready to display it
  _CLWindow = new CLGLWindow(800, 600,//height, width
			     0, 0,//initPosition (x,y)
			     "Visualizer : " + nName,
			     tickFreq,
			     true);
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    static_cast<CLGLWindow&>(*_CLWindow).addRenderObj(spec->getCoilRenderObj());

  CoilMaster::getInstance().addWindow(_CLWindow);

  {
    const magnet::thread::ScopedLock lock(static_cast<CLGLWindow&>(*_CLWindow).getDestroyLock());
    if (!_CLWindow->isReady()) return;

    BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
      {
	spec->updateRenderObj(static_cast<CLGLWindow&>(*_CLWindow).getCLState());
      }
    
    std::ostringstream os;
    os << "t:" << Sim->dSysTime;
    _CLWindow.as<CLGLWindow>().setSimStatus1(os.str());
    os.str("");
    os << "Events:" << Sim->eventCount;
    _CLWindow.as<CLGLWindow>().setSimStatus2(os.str());
  }

  I_cout() << "Visualizer initialised\nOpenCL Plaftorm:" 
	  << static_cast<CLGLWindow&>(*_CLWindow).getCLState().getPlatform().getInfo<CL_PLATFORM_NAME>()
	  << "\nOpenCL Device:" 
	  << static_cast<CLGLWindow&>(*_CLWindow).getCLState().getDevice().getInfo<CL_DEVICE_NAME>();
}

void
SVisualizer::runEvent() const
{
  _updateTime = _CLWindow.as<CLGLWindow>().getUpdateInterval();
  
  double locdt = dt;
  dt += _updateTime;

  //Update test
  if (_CLWindow.as<CLGLWindow>().simupdateTick())
    {
      //Actually move forward the system time
      Sim->dSysTime += locdt;
      Sim->ptrScheduler->stream(locdt);      
      //dynamics must be updated first
      Sim->dynamics.stream(locdt);
      locdt += Sim->freestreamAcc;
      Sim->freestreamAcc = 0;

      if (_CLWindow.as<CLGLWindow>().dynamoParticleSync())
	Sim->dynamics.getLiouvillean().updateAllParticles();
      
      BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
	Ptr->eventUpdate(*this, NEventData(), locdt);

      {
	const magnet::thread::ScopedLock lock(static_cast<CLGLWindow&>(*_CLWindow).getDestroyLock());
	if (!_CLWindow.as<CLGLWindow>().isReady()) return;
	
	BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
	  spec->updateRenderObj(_CLWindow.as<CLGLWindow>().getCLState());

	_CLWindow.as<CLGLWindow>().flagNewData();
      }

      std::ostringstream os;
      os << "t:" << Sim->dSysTime;
      
      _CLWindow.as<CLGLWindow>().setSimStatus1(os.str());
      os.str("");
      os << "Events:" << Sim->eventCount;
      _CLWindow.as<CLGLWindow>().setSimStatus2(os.str());
    }
}

void 
SVisualizer::initialise(size_t nID)
{ ID = nID; }

#endif

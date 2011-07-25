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
#include <coil/RenderObj/Function.hpp>
#include <coil/RenderObj/Spheres.hpp>
#include <coil/RenderObj/Lines.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

SVisualizer::SVisualizer(dynamo::SimData* nSim, std::string nName, double tickFreq):
  System(nSim)
{
  _updateTime = tickFreq * Sim->dynamics.units().unitTime();
  dt = -HUGE_VAL;//We want to ensure we get at least one update before
		 //anything occurs in the system
  
  sysName = "Visualizer";

  //Build a window, ready to display it
  _CLWindow = new coil::CLGLWindow("Visualizer : " + nName, tickFreq, true);
  
  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    if (spec.typeTest<CoilRenderObj>())
      _CLWindow.as<coil::CLGLWindow>().addRenderObj(dynamic_cast<const CoilRenderObj&>(*spec).getCoilRenderObj());

  BOOST_FOREACH(magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
    {
      CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));

      if (obj != NULL)
	_CLWindow.as<coil::CLGLWindow>().addRenderObj(obj->getCoilRenderObj());
    }

  _coil.getInstance().addWindow(_CLWindow);

  {
    const magnet::thread::ScopedLock lock(_CLWindow.as<coil::CLGLWindow>().getDestroyLock());
    if (!_CLWindow->isReady()) return;
    
    BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
      dynamic_cast<const CoilRenderObj&>(*spec).updateRenderData(_CLWindow.as<coil::CLGLWindow>().getGLContext());
    
    BOOST_FOREACH(magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
      {
	CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));
	
	if (obj != NULL) obj->updateRenderData(_CLWindow.as<coil::CLGLWindow>().getGLContext());
      }
    
    std::ostringstream os;
    os << "t:" << Sim->dSysTime;
    _CLWindow.as<coil::CLGLWindow>().setSimStatus1(os.str());
    os.str("");
    os << "Events:" << Sim->eventCount;
    _CLWindow.as<coil::CLGLWindow>().setSimStatus2(os.str());
  }
  
  _lastUpdate = boost::posix_time::microsec_clock::local_time();

  dout << "Visualizer initialised\nOpenCL Plaftorm:" 
       << _CLWindow.as<coil::CLGLWindow>().getGLContext().getCLPlatform().getInfo<CL_PLATFORM_NAME>()
       << "\nOpenCL Device:" 
       << _CLWindow.as<coil::CLGLWindow>().getGLContext().getCLDevice().getInfo<CL_DEVICE_NAME>() << std::endl;
}

void
SVisualizer::runEvent() const
{
  _updateTime = _CLWindow.as<coil::CLGLWindow>().getUpdateInterval();
  if (dt == -HUGE_VAL) dt = 0;
  
  double locdt = dt;
  dt += _updateTime;
  
  //Update test
  if (_CLWindow.as<coil::CLGLWindow>().simupdateTick())
    {
      //Actually move forward the system time
      Sim->dSysTime += locdt;
      Sim->ptrScheduler->stream(locdt);      
      //dynamics must be updated first
      Sim->dynamics.stream(locdt);
      locdt += Sim->freestreamAcc;
      Sim->freestreamAcc = 0;
      
      if (_CLWindow.as<coil::CLGLWindow>().dynamoParticleSync())
	Sim->dynamics.getLiouvillean().updateAllParticles();
      
      BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
	Ptr->eventUpdate(*this, NEventData(), locdt);
      
      {
	const magnet::thread::ScopedLock lock(_CLWindow.as<coil::CLGLWindow>().getDestroyLock());
	if (!_CLWindow.as<coil::CLGLWindow>().isReady()) return;
	
	BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
	  dynamic_cast<const CoilRenderObj&>(*spec).updateRenderData(_CLWindow.as<coil::CLGLWindow>().getGLContext());
	
	BOOST_FOREACH(magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
	  {
	    CoilRenderObj* obj = dynamic_cast<CoilRenderObj*>(&(*local));
	    
	    if (obj != NULL) obj->updateRenderData(_CLWindow.as<coil::CLGLWindow>().getGLContext());
	  }

	_CLWindow.as<coil::CLGLWindow>().flagNewData();
      }
      std::ostringstream os;
      os << "t:" << Sim->dSysTime;
      
      _CLWindow.as<coil::CLGLWindow>().setSimStatus1(os.str());
      os.str("");
      os << "Events:" << Sim->eventCount;
      _CLWindow.as<coil::CLGLWindow>().setSimStatus2(os.str());
    }

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

#endif

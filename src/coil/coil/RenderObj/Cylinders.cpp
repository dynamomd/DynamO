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
#include <coil/coilMaster.hpp>
#include "Cylinders.hpp"

namespace coil {
  RCylinders::RCylinders(size_t N, std::string name):
    RenderObj(name),
    _N(N)
  {}

  RCylinders::~RCylinders() { deinit(); }

  void 
  RCylinders::init(const magnet::thread::RefPtr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue);
    Cylinders::init(_N);
    _origin.resize(3 * _N);
    _orientation.resize(4 * _N);
    _scale.resize(3 * _N);
  }

  void 
  RCylinders::glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam)
  { Cylinders::glRender(); }

  void 
  RCylinders::deinit()
  { Cylinders::deinit(); }

  void 
  RCylinders::notifyDataUpdate()
  {
    CoilRegister::getCoilInstance().getTaskQueue()
      .queueTask(magnet::function::Task::makeTask(&RCylinders::dataUpdateWorker,this));
  }

  void 
  RCylinders::dataUpdateWorker()
  {
    _positionData = _origin;
    _orientationData = _orientation;
    _scalingData = _scale;
  }
}

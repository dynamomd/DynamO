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
#pragma once
#include "RenderObj.hpp"
#include <magnet/GL/objects/cylinders.hpp>
#include <vector>

namespace coil {
  class RCylinders : public RenderObj, public magnet::GL::objects::Cylinders
  {
  public:
    RCylinders(size_t N, std::string name);
    ~RCylinders();

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam);

    virtual void init(const magnet::thread::RefPtr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();

    std::vector<GLfloat>& getOriginData() { return _origin; }
    std::vector<GLfloat>& getOrientationData() { return _orientation; }
    std::vector<GLfloat>& getScaleData() { return _scale; }
    
    void notifyDataUpdate();

  protected:    
    void dataUpdateWorker();

    size_t _N;
    std::vector<GLfloat> _origin;
    std::vector<GLfloat> _orientation;
    std::vector<GLfloat> _scale;
  };
}

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

    virtual void glRender();
    virtual void initOpenGL();

    virtual void releaseCLGLResources();
  protected:
    size_t _N;
  };
}

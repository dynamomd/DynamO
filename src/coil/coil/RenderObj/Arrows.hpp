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
#pragma once

#include <coil/RenderObj/Lines.hpp>
#include <magnet/thread/mutex.hpp>

class RArrows : public RLines
{
public:
  RArrows(size_t N, std::string name);

  virtual void initOpenCL();
  virtual void initOpenGL();
  virtual void clTick(const magnet::GL::viewPort&);

  cl::Buffer& getPointData() { return _pointData; }
  cl::Buffer& getDirectionData() { return _directionData; }

protected:
  cl::Buffer _pointData;
  cl::Buffer _directionData;
  cl::Program _program;
  cl::Kernel _kernel;
  cl::KernelFunctor _kernelFunc;
};

/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "../clWindow.hpp"
#include "Triangles.hpp"
#include <time.h>


class RTTestWaves : public RTriangles
{
public:
  RTTestWaves(size_t N, float Yoffset);

  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&);

  void initOpenGL();
  void initOpenCL(magnet::CL::CLGLState&);

protected:
  static const std::string kernelsrc;

  cl::Kernel kernel;
  timespec startTime;

  size_t _N;
  float _Yoffset;

  cl::GLBuffer _clbuf_Normals;
};

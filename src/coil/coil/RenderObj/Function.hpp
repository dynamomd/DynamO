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


class RFunction : public RTriangles
{
public:
  RFunction(size_t N, Vector origin, Vector axis1,
	    Vector axis2, Vector axis3,
	    cl_float functionOriginX, 
	    cl_float functionOriginY,
	    cl_float functionRangeX,
	    cl_float functionRangeY,
	    bool drawAxis,
	    bool staticShape);

  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&);

  void initOpenGL();
  void initOpenCL(magnet::CL::CLGLState&);

  virtual void glRender();

  inline void setDrawAxis(bool val) { _drawAxis = val; }
  inline void setStaticShape(bool val) { _staticShape = val; }

protected:
  static const std::string kernelsrc;

  cl::Kernel _kernel;
  cl::KernelFunctor _kernelFunc;

  cl::Program _program;
  timespec startTime;

  cl_uint _N;

  Vector _origin;
  Vector _axis1;
  Vector _axis2;     
  Vector _axis3;

  cl_float3 _cl_origin;
  cl_float3 _cl_axis1;
  cl_float3 _cl_axis2;
  cl_float3 _cl_axis3;

  cl_float2 _functionOrigin;
  cl_float2 _functionRange;

  bool _drawAxis;
  bool _staticShape;
};

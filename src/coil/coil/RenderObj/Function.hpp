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

#include "Triangles.hpp"
#include <time.h>
#include <magnet/math/vector.hpp>

namespace coil {
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
	      bool staticShape,
	      std::string name,
	      std::string function = "f = pos.x * pos.y * native_sin(t);\n",
	      std::string normalCalc = "normal = normalize((float4)(pos.y * native_sin(t), pos.x * native_sin(t),1,0));\n",
	      std::string colorCalc = "\n");

    virtual void clTick(magnet::GL::Camera& cam) { clTick(); }

    void initOpenGL();
    void initOpenCL();

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam);

    inline void setDrawAxis(bool val) { _drawAxis = val; }
    inline void setStaticShape(bool val) { _staticShape = val; }

    virtual void initPicking(cl_uint& offset);
    virtual void pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam);
    virtual void finishPicking(cl_uint& offset, const cl_uint val);

    void setConstantA(cl_float val) { _A = val; }

    virtual void clTick(const magnet::GL::Camera&) {}

  protected:
    void clTick();
    std::string genKernelSrc();

    cl::Kernel _kernel;
    cl::Kernel _pickKernel;
    cl::KernelFunctor _kernelFunc;
    cl::KernelFunctor _pickFunc;
    std::string _kernelsrc;

    cl::Program _program;
    timespec startTime;

    cl_uint _N;
    cl_float tempo;

    Vector _origin;
    Vector _axis1;
    Vector _axis2;     
    Vector _axis3;

    cl_float4 _cl_origin;
    cl_float4 _cl_axis1;
    cl_float4 _cl_axis2;
    cl_float4 _cl_axis3;

    cl_float2 _functionOrigin;
    cl_float2 _functionRange;

    bool _drawAxis;
    bool _staticShape;
  
    std::string _function;
    std::string _normalCalc;
    std::string _colorCalc;

    volatile cl_float _A;
  };
}

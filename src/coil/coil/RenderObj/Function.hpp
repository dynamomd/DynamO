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
#include <magnet/CL/detail/program.hpp>

#define STRINGIFY(A) #A


namespace coil {
  class RFunction : public RTriangles
  {
    /* \brief An OpenCL program which converts a function into a triangle mesh.
     */
    class PlotProgram: public magnet::CL::detail::Program
    {
    public:
      PlotProgram(std::string function,
		  std::string normalCalc,
		  std::string colorCalc):
	_function(function),
	_normalCalc(normalCalc),
	_colorCalc(colorCalc)    
      {}

      virtual std::string initKernelSrc()
      {
	return std::string(STRINGIFY(
__kernel void
FunctionRenderKernel(__global float * positions,
		     __global uchar4 * colors,
		     __global float * normals,
		     float t,
		     float2 functionOrigin,
		     float2 functionRange,
		     float4 axis1,
		     float4 axis2,
		     float4 axis3,
		     float4 origin,
		     uint N, float A)
{
  positions += 3 * get_global_id(0);
  normals += 3 * get_global_id(0);
  colors += get_global_id(0);

  float2 normPos = (float2)(get_global_id(0) % N, get_global_id(0) / N);
  normPos /= N;

  float2 pos = normPos * functionRange + functionOrigin;

  float f; 
  )) + _function + std::string(STRINGIFY(
  float4 vertexPosition = normPos.x * axis1 + normPos.y * axis2 + f * axis3 + origin;
  
  positions[0] = vertexPosition.x;
  positions[1] = vertexPosition.y;
  positions[2] = vertexPosition.z;

  float4 normal;
  )) + _normalCalc + std::string(STRINGIFY(
  normal *= (float4)(functionRange * length(axis3) , 1.0f / length(axis3), 0);
  
  float4 rotatedNormal 
  = normalize(normal.x * axis1 +
	      normal.y * axis2 +
	      normal.z * axis3
	      );

  normals[0] = rotatedNormal.x;
  normals[1] = rotatedNormal.y;
  normals[2] = rotatedNormal.z;
  )) + _colorCalc + std::string(STRINGIFY(
  }					
__kernel void
FunctionPickKernel(__global uint * colors, uint offset)
{
  colors[get_global_id(0)] = get_global_id(0) + offset;
}
  ));
      }
      
    protected:

      std::string _function;
      std::string _normalCalc;
      std::string _colorCalc;
    };

  public:
    RFunction(std::string name,
	      size_t N = 10, Vector origin = Vector(-10,-1.0,-10), Vector axis1 = Vector(20,0,0),
	      Vector axis2 = Vector(0,0,20), Vector axis3 = Vector(0,1,0),
	      cl_float functionOriginX = 1.0,
	      cl_float functionOriginY = 1.0,
	      cl_float functionRangeX = 1.0,
	      cl_float functionRangeY = 1.0,
	      bool drawAxis = false,
	      bool staticShape = true,
	      std::string function = "f = 0.0 /*pos.x * pos.y * native_sin(t)*/;\n",
	      std::string normalCalc = "normal = normalize((float4)(0.0, 0.0, 1.0, 0.0));\n",
	      std::string colorCalc = "\n");

    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    virtual void glRender(const magnet::GL::Camera& cam, RenderMode mode = DEFAULT);

    inline void setDrawAxis(bool val) { _drawAxis = val; }
    inline void setStaticShape(bool val) { _staticShape = val; }

    // Cannot do picking for functions, as the triangle mesh will interpolate colors and screw it all up!
//    virtual void pickingRender(const magnet::GL::Camera& cam, uint32_t& offset);
//    virtual bool finishPicking(uint32_t& offset, const uint32_t val);

    void setConstantA(cl_float val) { _A = val; }

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    virtual bool deletable() { return false; }

  protected:
    void clTick();

    cl::Kernel _kernel;
    cl::Kernel _pickKernel;
    cl::KernelFunctor _kernelFunc;
    cl::KernelFunctor _pickFunc;
    std::string _kernelsrc;

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
  
    PlotProgram _program;
      
    volatile cl_float _A;
  };
}

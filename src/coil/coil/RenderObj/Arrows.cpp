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
#include "Arrows.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

#define STRINGIFY(A) #A

static cl_float4 getclVec(Vector vec)
{ 
  cl_float4 clvec;
  clvec.x = vec[0];
  clvec.y = vec[1];
  clvec.z = vec[2];
  clvec.w = 0;
  return clvec;
}

static const std::string lineKernelSource = STRINGIFY(
__kernel void
LineRenderKernel(const __global float* pointData,
		 const __global float* directionData,
		 __global float * vertexBuffer,
		 float4 camPos, float4 camdir, float4 camup,
		 uint Nlines)
{
  pointData += get_global_id(0) * 3;
  directionData += get_global_id(0) * 3;

//  float3 pos = *((const __global float3*)pointData);
//  float3 axis = *((const __global float3*)lineData);
}
						      );
RArrows::RArrows(size_t N):
  RLines(3*N) //3 lines per arrow
{}

void 
RArrows::initOpenCL(magnet::CL::CLGLState& CLState)
{
  try {
    _lineDataLock.lock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially lock the line data.";
    }

  RLines::initOpenCL(CLState);
  
  //Build buffer for line data
  _pointData = cl::Buffer(CLState.getContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY,
			  sizeof(cl_float) *  _N);

  _directionData = cl::Buffer(CLState.getContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY,
			      sizeof(cl_float) *  _N);

  //Build render kernel
  cl::Program::Sources kernelSource;
  kernelSource.push_back(std::pair<const char*, ::size_t>
			 (lineKernelSource.c_str(), lineKernelSource.size()));
  
  _program = cl::Program(CLState.getContext(), kernelSource);

  try {
    _program.build(std::vector<cl::Device>(1, CLState.getDevice()));
  } catch(cl::Error& err) {
    
    std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLState.getDevice());
    
    std::cout << "Compilation failed for device " <<
      CLState.getDevice().getInfo<CL_DEVICE_NAME>()
	      << "\nBuild Log:" << msg;
    
    throw;
  }
  
  _kernel = cl::Kernel(_program, "LineRenderKernel");

  cl_uint paddedN = (((_N/3) + 255)/256) * 256;

  _kernelFunc = _kernel.bind(CLState.getCommandQueue(), cl::NDRange(paddedN), cl::NDRange(256));
    
  try {
    _lineDataLock.unlock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially (UN)lock the line data.";
    }
}

void 
RArrows::clTick(magnet::CL::CLGLState& CLState, const magnet::GL::viewPort& _viewPortInfo)
{
  cl_float4 campos = getclVec(_viewPortInfo._position);
  cl_float4 camdir = getclVec(_viewPortInfo._cameraDirection);
  cl_float4 camup = getclVec(_viewPortInfo._cameraUp);

  //Aqquire GL buffer objects
  _clbuf_Positions.acquire(CLState.getCommandQueue());
  
  cl_uint NArrows = _N / 3;

  //Generate the sort data
  _kernelFunc(_pointData, _directionData, (cl::Buffer)_clbuf_Positions, 
	      campos, camdir, camup, NArrows);

  //Release resources
  _clbuf_Positions.release(CLState.getCommandQueue());
}

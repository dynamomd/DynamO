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
#include "Arrows.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>
#include <magnet/HSV.hpp>
#include <sstream>

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
		 float4 camPos, uint Nlines)
{
  //Position data
  if (get_global_id(0) >= Nlines) return;
  
  pointData += get_global_id(0) * 3;
  directionData += get_global_id(0) * 3;
  
  vertexBuffer += 4 * 3 * get_global_id(0); 
  
  float3 pos ;
  pos.x = pointData[0];
  pos.y = pointData[1];
  pos.z = pointData[2];
  
  float3 dir ;
  dir.x = directionData[0];
  dir.y = directionData[1];
  dir.z = directionData[2];
  
  float3 point = pos - 0.5f * dir;
  
  //Arrow Bottom
  vertexBuffer[0] = point.x;
  vertexBuffer[1] = point.y;
  vertexBuffer[2] = point.z;
  
  //Arrow Head
  point = pos + 0.5f * dir;
  vertexBuffer[3] = point.x;
  vertexBuffer[4] = point.y;
  vertexBuffer[5] = point.z;

  float3 pointToView = point - camPos.xyz;
  float3 sidesVec = normalize(cross(pointToView, dir));
  
  //Arrow verts
  point = pos + 0.3f * dir + 0.1 * length(dir) * sidesVec;
  vertexBuffer[6] = point.x;
  vertexBuffer[7] = point.y;
  vertexBuffer[8] = point.z;

  point = pos + 0.3f * dir - 0.1 * length(dir) * sidesVec;
  vertexBuffer[9] = point.x;
  vertexBuffer[10] = point.y;
  vertexBuffer[11] = point.z;
}
						      );
RArrows::RArrows(size_t N, std::string name):
  RLines(N, name)
{}

void 
RArrows::initOpenGL()
{
  {//Setup initial vertex positions, arrows have 4 verts
    std::vector<float> VertexPos(3 * _N * 4, 0.0);
    for (size_t i(0); i < _N; ++i)
      { 
	//Base
	VertexPos[12*i+0] = i * 1.0f / _N;
	VertexPos[12*i+1] = i * 1.0f / _N;
	VertexPos[12*i+2] = i * 1.0f / _N;

	//Head
	VertexPos[12*i+3] = i * 1.0f / _N;
	VertexPos[12*i+4] = (i + 0.5f) * 1.0f / _N;
	VertexPos[12*i+5] = i * 1.0f / _N;

	//Side 1
	VertexPos[12*i+6] = (i + 0.10f) * 1.0f / _N;
	VertexPos[12*i+7] = (i + 0.35f) * 1.0f / _N;
	VertexPos[12*i+8] = i * 1.0f / _N;
	
	//Side 2
	VertexPos[12*i+ 9] = (i - 0.10f) * 1.0f / _N;
	VertexPos[12*i+10] = (i + 0.35f) * 1.0f / _N;
	VertexPos[12*i+11] = i * 1.0f / _N;
      }
    setGLPositions(VertexPos);
  }
  
  {//4 vertexes per line
    std::vector<cl_uchar4> VertexColor(4 * _N);
    
    for (size_t icol = 0; icol < _N; ++icol)
      for (size_t jcol = 0; jcol < 4; ++jcol)
	magnet::color::HSVtoRGB(VertexColor[4*icol+jcol],
				float(icol)/ _N);

    setGLColors(VertexColor);
  }

  {//Setup initial element data
    //3 line segments, 2 vertices each
    std::vector<int> ElementData(6 * _N, 0);

    for (size_t i(0); i < _N; ++i)
      {
	//Base-head
	ElementData[6 * i + 0] = 4 * i + 0;
	ElementData[6 * i + 1] = 4 * i + 1;

	//head-side1
	ElementData[6 * i + 2] = 4 * i + 1;
	ElementData[6 * i + 3] = 4 * i + 2;

	//head-side2
	ElementData[6 * i + 4] = 4 * i + 1;
	ElementData[6 * i + 5] = 4 * i + 3;
      }
    
    setGLElements(ElementData);
  }
  
}

void 
RArrows::initOpenCL()
{
  RLines::initOpenCL();
  
  //Build buffer for line data
  _pointData = cl::Buffer(_CLState->getContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY,
			  sizeof(cl_float) *  _N * 3);

  _directionData = cl::Buffer(_CLState->getContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY,
			      sizeof(cl_float) *  _N * 3);

  //Build render kernel
  std::stringstream fullSource;

  fullSource << magnet::color::getOpenCLHSV();
  fullSource << lineKernelSource;
  
  //Need to make the c_str() point to a valid data area, so copy the string
  std::string finalSource = fullSource.str();

  cl::Program::Sources kernelSource;
  kernelSource.push_back(std::pair<const char*, ::size_t>
			 (finalSource.c_str(), finalSource.size()));

  _program = cl::Program(_CLState->getContext(), kernelSource);

  try {
    _program.build(std::vector<cl::Device>(1, _CLState->getDevice()));
  } catch(cl::Error& err) {
    
    std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_CLState->getDevice());
    
    std::cout << "Compilation failed for device " <<
      _CLState->getDevice().getInfo<CL_DEVICE_NAME>()
	      << "\nBuild Log:" << msg;
    
    throw;
  }
  
  _kernel = cl::Kernel(_program, "LineRenderKernel");

  cl_uint paddedN = ((_N + 255) / 256) * 256;

  _kernelFunc = _kernel.bind(_CLState->getCommandQueue(), cl::NDRange(paddedN), cl::NDRange(256));
}


void 
RArrows::clTick(const magnet::GL::viewPort& _viewPortInfo)
{
  cl_float4 campos = getclVec(_viewPortInfo._position);

  //Aqquire GL buffer objects
  _clbuf_Positions.acquire(_CLState->getCommandQueue());
  
  cl_uint NArrows = _N;

  //Generate the sort data
  _kernelFunc(_pointData, _directionData, (cl::Buffer)_clbuf_Positions, 
	      campos, NArrows);
  
  //Release resources
  _clbuf_Positions.release(_CLState->getCommandQueue());
}

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

#include "Function.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

#define STRINGIFY(A) #A

const std::string 
RFunction::kernelsrc = STRINGIFY(
__constant float decayrate = 2.5f;
__constant float invWaveLength = 40.0f;
__constant float freq = -4;

float func(float2 pos, float t)
{
  float r = native_sqrt(pos.x * pos.x + pos.y * pos.y);
  return native_exp( - decayrate * r) * native_sin(invWaveLength * r + freq * t);
}

float3 funcNormal(float2 pos, float t)
{
  float r = native_sqrt(pos.x * pos.x + pos.y * pos.y);

  float dfodr = native_exp(- decayrate * r) 
    * (invWaveLength * native_cos(r * invWaveLength + freq * t)
       + decayrate * native_sin(r * invWaveLength + freq * t));

  return normalize((float3)(-dfodr * pos.x / r, -dfodr * pos.y / r,1));
}

__kernel void
TestWaveKernel(__global float * positions,
	       __global float * colors,
	       __global float * normals,
	       float t,
	       float2 functionOrigin,
	       float2 functionRange,
	       float3 axis1,
	       float3 axis2,
	       float3 axis3,
	       float3 origin,
	       uint N
	       )
{
  float2 normPos = (float2)(get_global_id(0) % N, get_global_id(0) / N);
  normPos /= N;

  float2 functionPosition = normPos * functionRange + functionOrigin;

  float val = func(functionPosition, t);
  
  float3 vertexPosition = normPos.x * axis1 + normPos.y * axis2 + val * axis3 + origin;

  const uint id = get_global_id(0);

  positions[3*id+0] = vertexPosition.x;
  positions[3*id+1] = vertexPosition.y;
  positions[3*id+2] = vertexPosition.z;

  //colors[id] = (float4)(clamp(val, 0.0f, 1.0f),0,0,1); 

  float3 normal = funcNormal(functionPosition, t);
  float3 rotatedNormal = normalize(normal.x * axis1 + normal.y * axis2 + normal.z * axis3);
  normals[3*id+0] = rotatedNormal.x;
  normals[3*id+1] = rotatedNormal.y;
  normals[3*id+2] = rotatedNormal.z;
}
);

RFunction::RFunction(size_t N, Vector origin, Vector axis1,
		     Vector axis2, Vector axis3, cl_float functionOriginX,
		     cl_float functionOriginY, cl_float functionRangeX, 
		     cl_float functionRangeY, bool drawAxis, bool staticShape):
  _N(N),
  _origin(origin),
  _axis1(axis1),
  _axis2(axis2),
  _axis3(axis3),
  _drawAxis(drawAxis),
  _staticShape(staticShape)
{
  //Copy to the cl types
  for (size_t i(0); i < 3; ++i)
    {
      _cl_origin.s[i] = origin[i];
      _cl_axis1.s[i] = axis1[i];
      _cl_axis2.s[i] = axis2[i];
      _cl_axis3.s[i] = axis3[i];
    }

  _functionOrigin.s[0] = functionOriginX;
  _functionOrigin.s[1] = functionOriginY;
  _functionRange.s[0] = functionRangeX;
  _functionRange.s[1] = functionRangeY;

  //Change N to a multiple of 16 so the workgroup size of 256 always fits
  _N = ((_N +15) / 16) * 16;
}

void 
RFunction::initOpenGL()
{}

void 
RFunction::initOpenCL(magnet::CL::CLGLState& CLState)
{
  {//Setup initial vertex positions
    float spacingFactor = 1.0 / (_N + 0.5f);

    std::vector<float> VertexPos(3 * _N * _N, 0.0);  
    for (size_t i = 0; i < _N; i++)
      {       
	for (size_t j = 0; j < _N; j++)
	  {
	    Vector pos = (i * _axis1 + j * _axis2) * spacingFactor + _origin;
	    
	    VertexPos[0 + 3 * (i + _N * j)] = pos[0];	
	    VertexPos[1 + 3 * (i + _N * j)] = pos[1];
	    VertexPos[2 + 3 * (i + _N * j)] = pos[2];
	  }
      }
    setGLPositions(VertexPos);
  }

  {//Setup inital normal vectors
    std::vector<float> VertexNormals(3 * _N * _N, 0.0);
    Vector normal = _axis1 ^ _axis2;	

    for (size_t i = 0; i < _N; i++)
      {
	for (size_t j = 0; j < _N; j++)
	  {
	    VertexNormals[0 + 3 * (i + _N * j)] = normal[0];
	    VertexNormals[1 + 3 * (i + _N * j)] = normal[1];
	    VertexNormals[2 + 3 * (i + _N * j)] = normal[2];
	  }
      }
    setGLNormals(VertexNormals);
  }

  {//Setup initial Colors
    std::vector<float> VertexColor(4 * _N * _N, 0.0);
    for (size_t i = 0; i < _N; i++)
      {       
	for (size_t j = 0; j < _N; j++)
	  {
	    VertexColor[0 + 4 * (i + _N * j)] = 1.0f;
	    VertexColor[1 + 4 * (i + _N * j)] = 1.0f;
	    VertexColor[2 + 4 * (i + _N * j)] = 1.0f;
	    VertexColor[3 + 4 * (i + _N * j)] = 1.0f;
	  }       
      }
    setGLColors(VertexColor);
  }
   
  {//Setup initial element data
    std::vector<int> ElementData(3*2*(_N-1)*(_N-1), 0.0);
    for  (size_t i = 0; i < _N - 1; i++)
      {
	for (size_t j = 0; j < _N - 1; j++)
	  {
	    ElementData[6 * (i + (_N - 1) * j) + 0] = i + _N * j;
	    ElementData[6 * (i + (_N - 1) * j) + 1] = i + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 2] = i + 1 + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 3] = i + _N * j;
	    ElementData[6 * (i + (_N - 1) * j) + 4] = i + 1 + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 5] = i + 1 + _N * j;
	  }
      }
    setGLElements(ElementData);
  }
  
  cl::Program::Sources kernelSource;
  kernelSource.push_back(std::pair<const char*, ::size_t>(kernelsrc.c_str(), kernelsrc.size()));
  
  _program = cl::Program(CLState.getCommandQueue().getInfo<CL_QUEUE_CONTEXT>(), kernelSource);
  
  try 
    {
      _program.build(std::vector<cl::Device>(1, CLState.getDevice()));
    }
  catch(cl::Error& err) 
    {    
      std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLState.getDevice());
      
      std::cout << "Compilation failed for device " <<
	CLState.getDevice().getInfo<CL_DEVICE_NAME>()
		<< "\nBuild Log:" << msg;
      throw;
    }
  
  _kernel = cl::Kernel(_program, "TestWaveKernel");
  
  const size_t workgroupSize = 256;

  //N is a multiple of 16, so a workgroup size of 256 is always good
  _kernelFunc = _kernel.bind(CLState.getCommandQueue(), cl::NDRange(_N * _N), 
			     cl::NDRange(workgroupSize));

  clock_gettime(CLOCK_MONOTONIC, &startTime);
  glFinish();

  initOCLVertexBuffer(CLState.getContext());
  initOCLNormBuffer(CLState.getContext());
  initOCLColorBuffer(CLState.getContext());

  //Now do the first clTick if the shape is static!
  if (_staticShape)
    {
      _staticShape = false;
      clTick(CLState, magnet::GL::viewPort());
      _staticShape = true;      
    }
}


void 
RFunction::clTick(magnet::CL::CLGLState& CLState, const magnet::GL::viewPort& _viewPortInfo)
{
  if (_staticShape) return;

  timespec currTime;
  clock_gettime(CLOCK_MONOTONIC, &currTime);
  
  float tempo = float(currTime.tv_sec) - float(startTime.tv_sec)
    + 1e-9 * (float(currTime.tv_nsec) - float(startTime.tv_nsec));

  //Aqquire buffer objects
  _clbuf_Colors.acquire(CLState.getCommandQueue());
  _clbuf_Positions.acquire(CLState.getCommandQueue());
  _clbuf_Normals.acquire(CLState.getCommandQueue());
      
  //Run Kernel
  _kernelFunc((cl::Buffer)_clbuf_Positions, 
	      (cl::Buffer)_clbuf_Colors, 
	      (cl::Buffer)_clbuf_Normals, 
	      tempo,
	      _functionOrigin,
	      _functionRange,
	      _cl_axis1,
	      _cl_axis2,
	      _cl_axis3,
	      _cl_origin,
	      _N);
  
  //Release resources
  _clbuf_Normals.release(CLState.getCommandQueue());
  _clbuf_Colors.release(CLState.getCommandQueue());
  _clbuf_Positions.release(CLState.getCommandQueue());
}

void 
RFunction::glRender()
{
  RTriangles::glRender();

  //Draw the axis of the plotted function
  if (_drawAxis)
    {
      coil::glprimatives::drawArrow(_origin, _origin + _axis1);
      coil::glprimatives::drawArrow(_origin, _origin + _axis2);
      coil::glprimatives::drawArrow(_origin, _origin + _axis3);
    }
}

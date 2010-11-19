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

#include "TestWaves.hpp"
#include <iostream>

#define STRINGIFY(A) #A

const std::string 
RTTestWaves::kernelsrc = STRINGIFY(
__constant float decayrate = 2.5f;
__constant float invWaveLength = 40.0f;
__constant float freq = -4;

float wavefunc(float x, float z, float t)
{
  float r = native_sqrt(x * x + z * z);
  return native_exp( - decayrate * r) * native_sin(invWaveLength * r + freq * t);
}

float3 waveNormal(float x, float z, float t)
{
  float r = native_sqrt(x * x + z * z);

  float dfodr = native_exp(- decayrate * r) 
    * (invWaveLength * native_cos(r * invWaveLength + freq * t)
       + decayrate * native_sin(r * invWaveLength + freq * t));

  return normalize((float3)(- dfodr * x / r, 1, - dfodr * z / r));
}

__kernel void
TestWaveKernel(__global float * positions,
	       __global float * colors,
	       //__global float * normals,
	       float t, float Yoffset)
{
  int i = get_global_id(0);
  
  float x = positions[3*i];
  float z = positions[3*i+2];
  
  float val =  wavefunc(x + 0.7f, z, t) + wavefunc(x - 0.7f, z, t) + Yoffset;
  positions[3*i+1] = val;

  //float3 normal = 0.5 * (waveNormal(x + 0.7f, z, t) + waveNormal(x - 0.7f, z, t));
  //normals[3*i+0] = normal.x;
  //normals[3*i+1] = normal.y;
  //normals[3*i+2] = normal.z;

  colors[4*i+0] = clamp(val, 0.0f, 1.0f); 
  
//  colors[4*i+0] = clamp(dot(normal, (float3)(0,0,1)), 0.0f, 1.0f);
//  colors[4*i+1] = clamp(dot(normal, (float3)(0,0,1)), 0.0f, 1.0f);
//  colors[4*i+2] = clamp(dot(normal, (float3)(0,0,1)), 0.0f, 1.0f);
}
);

RTTestWaves::RTTestWaves(size_t N, float Yoffset):
  _N(N),
  _Yoffset(Yoffset)
{
}

void 
RTTestWaves::initOpenGL() 
{}

void 
RTTestWaves::initOpenCL(magnet::CL::CLGLState& CLState)
{
  {//Setup initial vertex positions
    std::vector<float> VertexPos(3 * _N * _N, 0.0);  
    for (size_t i = 0; i < _N; i++)
      {       
	for (size_t j = 0; j < _N; j++)
	  {
	    VertexPos[0 + 3 * (i + _N * j)] = 4*((float)i / (float)_N-0.5f);	
	    VertexPos[1 + 3 * (i + _N * j)] = 0.0;
	    VertexPos[2 + 3 * (i + _N * j)] = 4*((float)j / (float)_N-0.5f);
	  }       
      }
    setGLPositions(VertexPos);
    initOCLVertexBuffer(CLState.getContext());
  }

  {//Setup inital normal vectors
    std::vector<float> VertexNormals(3 * _N * _N, 0.0);
    for (size_t i = 0; i < _N; i++)
      {       
	for (size_t j = 0; j < _N; j++)
	  {
	    VertexNormals[0 + 3 * (i + _N * j)] = 0.0f;
	    VertexNormals[1 + 3 * (i + _N * j)] = 1.0f;
	    VertexNormals[2 + 3 * (i + _N * j)] = 0.0f;
	  }       
      }
    setGLNormals(VertexNormals);
    initOCLNormBuffer(CLState.getContext());
  }

  {//Setup initial Colors
    std::vector<float> VertexColor(4 * _N * _N, 0.0);
    for (size_t i = 0; i < _N; i++)
      {       
	for (size_t j = 0; j < _N; j++)
	  {
	    VertexColor[0 + 4 * (i + _N * j)] = 0.0f;
	    VertexColor[1 + 4 * (i + _N * j)] = 0.0f;
	    VertexColor[2 + 4 * (i + _N * j)] = (float)i/(float)(_N-1);
	    VertexColor[3 + 4 * (i + _N * j)] = 1.0f;
	  }       
      }
    setGLColors(VertexColor);
    initOCLColorBuffer(CLState.getContext());
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
  
  program = cl::Program(CLState.getCommandQueue().getInfo<CL_QUEUE_CONTEXT>(), kernelSource);
  
  try {
    program.build(std::vector<cl::Device>(1, CLState.getDevice()));
  } catch(cl::Error& err) {
    
    std::string msg = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLState.getDevice());
    
    std::cout << "Compilation failed for device " <<
      CLState.getDevice().getInfo<CL_DEVICE_NAME>()
	      << "\nBuild Log:" << msg;
    throw;
  }
  
  kernel = cl::Kernel(program, "TestWaveKernel");
  
  clock_gettime(CLOCK_MONOTONIC, &startTime);
}


void 
RTTestWaves::clTick(magnet::CL::CLGLState& CLState, const magnet::GL::viewPort& _viewPortInfo)
{
  //RTriangles::clTick(CLState, _viewPortInfo);

  cl::KernelFunctor kernelFunc = kernel.bind(CLState.getCommandQueue(), cl::NDRange(_N * _N), cl::NDRange(200));
  timespec currTime;
  clock_gettime(CLOCK_MONOTONIC, &currTime);
  
  float tempo = float(currTime.tv_sec) - float(startTime.tv_sec)
    + 1e-9 * (float(currTime.tv_nsec) - float(startTime.tv_nsec));

  //Aqquire buffer objects

  _clbuf_Colors.acquire(CLState.getCommandQueue());
  _clbuf_Positions.acquire(CLState.getCommandQueue());
  try {
    //_clbuf_Normals.acquire(CLState.getCommandQueue());
      
    //Run Kernel
    kernelFunc((cl::Buffer)_clbuf_Positions, 
	       (cl::Buffer)_clbuf_Colors, 
	       //(cl::Buffer)_clbuf_Normals, 
	       tempo, _Yoffset);
      
      //Release resources
    //_clbuf_Normals.release(CLState.getCommandQueue());
  } catch (cl::Error& err)
    {
      std::cerr << "OpenCL error: " << err.what() << "(" << err.err() <<
	")" << std::endl;
    }
  
  _clbuf_Colors.release(CLState.getCommandQueue());
  _clbuf_Positions.release(CLState.getCommandQueue());
}

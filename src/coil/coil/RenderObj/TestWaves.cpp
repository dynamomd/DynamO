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
__kernel void
TestWaveKernel(__global float * positions,__global float * cores, float t, float Yoffset)
{
  int i = get_global_id(0);
  
  float x = positions[3*i]+0.7f;
  float y = positions[3*i+2];
  float r = native_sqrt(x*x+y*y);
  
  float valor = native_exp(- r * 2.5f)*native_sin(40*r-4*t);
  x -= 1.4f;
  r = native_sqrt(x*x+y*y);
  valor += native_exp(- r * 1.5f)*native_sin(40*r-4*t);
  
  positions[3*i+1] = valor + Yoffset;
  cores[4*i] = clamp(valor,0,1); 
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
RTTestWaves::initOpenCL(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers)
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
    initOCLVertexBuffer(Context, hostTransfers);
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
    initOCLColorBuffer(Context, hostTransfers);
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
  
  cl::Program program(CmdQ.getInfo<CL_QUEUE_CONTEXT>(), kernelSource);
  
  std::string buildOptions;
  
  cl::Device clDevice = CmdQ.getInfo<CL_QUEUE_DEVICE>();
  try {
    program.build(std::vector<cl::Device>(1, clDevice), buildOptions.c_str());
  } catch(cl::Error& err) {
    
    std::string msg = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(Device);
    
    std::cout << "Compilation failed for device " <<
      Device.getInfo<CL_DEVICE_NAME>()
	      << "\nBuild Log:" << msg;
    
    throw;
  }
  
  kernel = cl::Kernel(program, "TestWaveKernel");
  
  clock_gettime(CLOCK_MONOTONIC, &startTime);
}


void 
RTTestWaves::clTick(cl::CommandQueue& CmdQ, cl::Context& Context)
{
  cl::KernelFunctor kernelFunc = kernel.bind(CmdQ, cl::NDRange(_N * _N), cl::NDRange(200));
  timespec currTime;
  clock_gettime(CLOCK_MONOTONIC, &currTime);
  
  float tempo = float(currTime.tv_sec) - float(startTime.tv_sec)
    + 1e-9 * (float(currTime.tv_nsec) - float(startTime.tv_nsec));

  //Aqquire buffer objects
  _clbuf_Colors.acquire(CmdQ);
  _clbuf_Positions.acquire(CmdQ);
  
  //Run Kernel
  kernelFunc((cl::Buffer)_clbuf_Positions, (cl::Buffer)_clbuf_Colors, tempo, _Yoffset);
  
  //Release resources
  _clbuf_Colors.release(CmdQ);
  _clbuf_Positions.release(CmdQ);
}

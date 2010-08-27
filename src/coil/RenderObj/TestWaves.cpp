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

RTTestWaves::RTTestWaves(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, 
			 bool hostTransfers, size_t N, float Yoffset):
  RTriangles(hostTransfers),
  _N(N),
  _Yoffset(Yoffset)
{
  {//Setup initial vertex positions
    std::vector<float> VertexPos(3 * N * N, 0.0);  
    for (size_t i = 0; i < N; i++)
      {       
	for (size_t j = 0; j < N; j++)
	  {
	    VertexPos[0 + 3 * (i + N * j)] = 4*((float)i / (float)N-0.5f);	
	    VertexPos[1 + 3 * (i + N * j)] = 0.0;
	    VertexPos[2 + 3 * (i + N * j)] = 4*((float)j / (float)N-0.5f);
	  }       
      }
    setGLPositions(VertexPos);
    initOCLVertexBuffer(Context);
  }

  {//Setup inital normal vectors
    std::vector<float> VertexNormals(3 * N * N, 0.0);
    for (size_t i = 0; i < N; i++)
      {       
	for (size_t j = 0; j < N; j++)
	  {
	    VertexNormals[0 + 3 * (i + N * j)] = 0.0f;
	    VertexNormals[1 + 3 * (i + N * j)] = 1.0f;
	    VertexNormals[2 + 3 * (i + N * j)] = 0.0f;
	  }       
      }
    setGLNormals(VertexNormals);
  }

  {//Setup initial Colors
    std::vector<float> VertexColor(4 * N * N, 0.0);
    for (size_t i = 0; i < N; i++)
      {       
	for (size_t j = 0; j < N; j++)
	  {
	    VertexColor[0 + 4 * (i + N * j)] = 0.0f;
	    VertexColor[1 + 4 * (i + N * j)] = 0.0f;
	    VertexColor[2 + 4 * (i + N * j)] = (float)i/(float)(N-1);
	    VertexColor[3 + 4 * (i + N * j)] = 1.0f;
	  }       
      }
    setGLColors(VertexColor);
    initOCLColorBuffer(Context);
  }
   
  {//Setup initial element data
    std::vector<int> ElementData(3*2*(N-1)*(N-1), 0.0);
    for  (size_t i = 0; i < N - 1; i++)
      {
	for (size_t j = 0; j < N - 1; j++)
	  {
	    ElementData[6 * (i + (N - 1) * j) + 0] = i + N * j;
	    ElementData[6 * (i + (N - 1) * j) + 1] = i + N * (j + 1);
	    ElementData[6 * (i + (N - 1) * j) + 2] = i + 1 + N * (j + 1);
	    ElementData[6 * (i + (N - 1) * j) + 3] = i + N * j;
	    ElementData[6 * (i + (N - 1) * j) + 4] = i + 1 + N * (j + 1);
	    ElementData[6 * (i + (N - 1) * j) + 5] = i + 1 + N * j;
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

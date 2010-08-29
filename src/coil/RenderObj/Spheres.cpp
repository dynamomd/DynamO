#include "Spheres.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <fstream>

#include "Primatives/Sphere.hpp"
#include "Spheres.clh"

struct  SortDataType { cl_uint ID; cl_float dist;};


RTSpheres::RTSpheres(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers,
		     const float& cameraX, const float& cameraY, const float& cameraZ,
		     size_t N, 
		     Sphere::SphereType type1, 
		     size_t order1,
		     Sphere::SphereType type2, 
		     size_t order2,
		     size_t nSphere1):
  RTriangles(hostTransfers),
  _N(N),
  primSphere1(type1, order1),
  primSphere2(type2, order2),
  _workgroupsize(0),
  _globalsize(0),
  _nSpheres1(nSphere1),
  _cameraX(cameraX),
  _cameraY(cameraY),
  _cameraZ(cameraZ)
{
  {
    size_t Ncuberoot = (size_t)std::pow(_N, 1.0/3.0);
    
    //We add one on here to allow vector loads and stores to occur
    _spherePositions = cl::Buffer(Context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, 
				  sizeof(cl_float) * (3 * _N + 1));

    //Work out the details for the bitonic sort
    _sortData = cl::Buffer(Context, CL_MEM_READ_WRITE, sizeof(SortDataType) * _N);
    
    _numStages = 0;
    for (size_t temp = _N; temp > 1; temp >>= 1)
      ++_numStages;

    _powerOfTwo = 1 << _numStages;
    if (_powerOfTwo < _N) { _powerOfTwo <<= 1; ++_numStages; }


    cl_float* Pos = (cl_float*)CmdQ.enqueueMapBuffer(_spherePositions, true, 
						     CL_MAP_WRITE, 0, 
						     3 * _N * sizeof(cl_float));
    
    //Generates positions on a simple cubic lattice
    for (size_t partID(0); partID < _N; ++partID)
      {
	Pos[3 * partID + 0] = (partID % Ncuberoot);
	Pos[3 * partID + 1] = ((partID / Ncuberoot)  % Ncuberoot);
	Pos[3 * partID + 2] = (partID / (Ncuberoot * Ncuberoot));
      }

    //Start copying this data to the graphics card
    CmdQ.enqueueUnmapMemObject(_spherePositions, (void*)Pos);
  }

  {//Setup initial vertex positions
    size_t nVerticeComp = 3 * (primSphere1.n_vertices * nSphere1 + primSphere2.n_vertices * (_N - nSphere1));
    std::vector<float> VertexPos(nVerticeComp, 0.0);

    setGLPositions(VertexPos);

    initOCLVertexBuffer(Context);
  }
  
  {//Setup inital normal vectors
    size_t nNormals1 = primSphere1.n_vertices * nSphere1;
    size_t nNormals2 = primSphere2.n_vertices * (_N - nSphere1);
    std::vector<float> VertexNormals(3 * (nNormals1 + nNormals2), 0.0);

    for (size_t i = 0; i < nSphere1; ++i)
      for (int j = 0; j < 3 * primSphere1.n_vertices; ++j)
	VertexNormals[3 * primSphere1.n_vertices * i + j] = primSphere1.vertices[j];
    
    for (size_t i = 0; i < (_N - nSphere1); ++i)
      for (int j = 0; j < 3 * primSphere2.n_vertices; ++j)
	VertexNormals[3 * (primSphere2.n_vertices * i + nNormals1) + j]
	  = primSphere2.vertices[j];

    setGLNormals(VertexNormals);
  }


  {//Setup initial Colors
    size_t nColors1 = 4 * primSphere1.n_vertices * nSphere1;
    size_t nColors2 = 4 * primSphere2.n_vertices * (_N - nSphere1);
    std::vector<float> VertexColor(nColors1 + nColors2, 1.0);

    setGLColors(VertexColor);
  }
   
 
  {//Setup initial element data
    size_t nElem1 = 3 * primSphere1.n_faces * nSphere1;
    size_t nElem2 = 3 * primSphere2.n_faces * (_N - nSphere1);

    std::vector<int> ElementData(nElem1 + nElem2, 0.0);
    
    for (size_t i = 0; i < nSphere1; ++i)
      for (int j = 0; j < 3 * primSphere1.n_faces; ++j)
	ElementData[3 * primSphere1.n_faces * i + j] 
	  = i * primSphere1.n_vertices + primSphere1.faces[j];

    for (size_t i = 0; i < _N - nSphere1; ++i)
      for (int j = 0; j < 3 * primSphere2.n_faces; ++j)
	ElementData[nElem1 + 3 * primSphere2.n_faces * i + j] 
	  = primSphere1.n_vertices * nSphere1 + i * primSphere2.n_vertices + primSphere2.faces[j];

    setGLElements(ElementData);
  }
  

  std::stringstream fullSource;
  
  //It is ideal if the workgroup size divides by 3(coords), 64
  //(warp/wave) AND the number of vertices per particle (not so important)

  //An Icosahedron, of order 0 (12), fits exactly into
  //3x32x2=192=12x16
  _workgroupsize = 2*32*3;
  _globalsize = _workgroupsize * (std::min((_N +_workgroupsize-1) / _workgroupsize, 
					   _workgroupsize*(9216 / _workgroupsize)));

  fullSource << "#define WORKGROUP_SIZE " << _workgroupsize << "\n";
  
  fullSource << sphereKernelSource;
  
  //Need to make the c_str() point to a valid data area, so copy the string
  std::string finalSource = fullSource.str();

  cl::Program::Sources kernelSource;
  kernelSource.push_back(std::pair<const char*, ::size_t>
			 (finalSource.c_str(), finalSource.size()));
  
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
  
  _renderKernel = cl::Kernel(program, "SphereRenderKernel");
  _sortDataKernel = cl::Kernel(program, "GenerateData");
  _sortKernel = cl::Kernel(program, "sphereBitonicSort");

  _primativeVertices1 = cl::Buffer(Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   sizeof(cl_float) * 3 * primSphere1.n_vertices,
				   primSphere1.vertices);

  _primativeVertices2 = cl::Buffer(Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   sizeof(cl_float) * 3 * primSphere2.n_vertices,
				   primSphere2.vertices);

  clTick(CmdQ, Context);
}

void 
RTSpheres::clTick(cl::CommandQueue& CmdQ, cl::Context& Context)
{
  cl::KernelFunctor renderKernelFunc = _renderKernel.bind(CmdQ, cl::NDRange(_globalsize), cl::NDRange(_workgroupsize));
  cl::KernelFunctor sortDataKernelFunc = _sortDataKernel.bind(CmdQ, cl::NDRange(_globalsize), cl::NDRange(_workgroupsize));
  cl::KernelFunctor sortKernelFunc = _sortKernel.bind(CmdQ, cl::NDRange(_powerOfTwo), cl::NDRange(256));

  //Generate the sort data
  sortDataKernelFunc(_spherePositions, _sortData, 
		     (cl_float)_cameraX, (cl_float)_cameraY, 
		     (cl_float)_cameraZ, _N);

  //Now sort the data
  for (cl_uint stage = 0; stage < _numStages-1; ++stage)
    for (cl_uint stagePass = 0; stagePass < stage+1; ++stagePass)
      sortKernelFunc(_sortData, stage, stagePass, _N, 0);

  {
    cl_uint stage = _numStages-1;	
    for (cl_uint stagePass = 0; stagePass < _numStages; ++stagePass)
      sortKernelFunc(_sortData, stage, stagePass, _N, 1);
  }

  //Aqquire GL buffer objects
  _clbuf_Positions.acquire(CmdQ);

  //Finally, run rendering kernel
  renderKernelFunc(_spherePositions, (cl::Buffer)_clbuf_Positions, _primativeVertices1, 
		   primSphere1.n_vertices, 0, _nSpheres1, 0, _sortData);

  renderKernelFunc(_spherePositions, (cl::Buffer)_clbuf_Positions, _primativeVertices2, 
		   primSphere2.n_vertices, _nSpheres1, _N, 
		   3 * _nSpheres1 * (primSphere1.n_vertices - primSphere2.n_vertices),
		   _sortData);

  //Release resources
  _clbuf_Positions.release(CmdQ);
}

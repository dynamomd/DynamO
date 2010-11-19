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

#include "Spheres.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <fstream>

#include "Spheres.clh"
#include <errno.h>

struct  SortDataType { cl_uint ID; cl_float dist; };

cl_float4 getclVec(Vector vec)
{ 
  cl_float4 clvec;
  clvec.x = vec[0];
  clvec.y = vec[1];
  clvec.z = vec[2];
  clvec.w = 0;
  return clvec;
}


RTSpheres::RTSpheres(size_t N):
  _N(N),
  _frameCount(0),
  _sortFrequency(1),
  _workgroupsize(0),
  _globalsize(0)
{
  //Work computer test render
  size_t spheres_rendered = 0;

  size_t stage_spheres = std::min(size_t(10ul), N - spheres_rendered);
  if (stage_spheres)
    {
      _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 2, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(size_t(1000ul), N - spheres_rendered);
  if (stage_spheres)
    {
      _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 1, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(size_t(10000ul), N - spheres_rendered);
  if (stage_spheres)
    {
      _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(size_t(200000ul), N - spheres_rendered);
  if (stage_spheres)
    {
      _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::octahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = N - spheres_rendered;
  if (stage_spheres)
    _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::tetrahedron, 0, N - spheres_rendered));
}

RTSpheres::~RTSpheres()
{
  //We get a lock to ensure we don't destroy during a data transfer
  magnet::thread::ScopedLock tmp(_sphereDataLock);
}

void
RTSpheres::initOpenCL(magnet::CL::CLGLState& CLState)
{
  try {
    _sphereDataLock.lock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially lock the sphere data.";
    }

  {
    _spherePositions = cl::Buffer(CLState.getContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, 
				  sizeof(cl_float4) *  _N);

    //We must pad the sort data out to a multiple of 1024

    size_t paddedN = ((_N + 1023)/1024) * 1024;
    
    _sortKeys = cl::Buffer(CLState.getContext(), CL_MEM_READ_WRITE, sizeof(float) * paddedN);
    _sortData = cl::Buffer(CLState.getContext(), CL_MEM_READ_WRITE, sizeof(cl_uint) * paddedN);
    
    cl_float4* Pos = (cl_float4*)CLState.getCommandQueue().enqueueMapBuffer(_spherePositions, true, 
									    CL_MAP_WRITE, 0, 
									    _N * sizeof(cl_float4));

    const float density = 0.1;

    cl_float particleDiam = std::pow(1 * density / _N, float(1.0 / 3.0));
    
    //Generates positions on a simple cubic lattice
    for (size_t partID(0); partID < _N; ++partID)
      {
	Pos[partID].x = ((1.0 * rand()) / RAND_MAX) - 0.5;
	Pos[partID].y = ((1.0 * rand()) / RAND_MAX) - 0.5;
	Pos[partID].z = ((1.0 * rand()) / RAND_MAX) - 0.5;

	Pos[partID].w = particleDiam * 0.5;
      }

    //Start copying this data to the graphics card
    CLState.getCommandQueue().enqueueUnmapMemObject(_spherePositions, (void*)Pos);
  }

  {//Setup initial vertex positions
    size_t nVertice = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nVertice += iPtr->_type.getVertexCount() * iPtr->_nSpheres;

    std::vector<float> VertexPos(3 * nVertice, 0.0);
    setGLPositions(VertexPos);
    initOCLVertexBuffer(CLState.getContext());
  }
  
  {//Setup inital normal vectors
    size_t nNormals = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nNormals += iPtr->_type.getVertexCount() * iPtr->_nSpheres;

    std::vector<float> VertexNormals(3 * nNormals, 0.0);

    nNormals = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	for (size_t i = 0; i < iPtr->_nSpheres; ++i)
	  for (int j = 0; j < 3 * iPtr->_type.getVertexCount(); ++j)
	    VertexNormals[nNormals + 3 * iPtr->_type.getVertexCount() * i + j] = iPtr->_type.getVertices()[j];

	nNormals += 3 * iPtr->_nSpheres * iPtr->_type.getVertexCount();
      }
    
    setGLNormals(VertexNormals);
  }


  {//Setup initial Colors
    size_t nColors = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nColors += iPtr->_type.getVertexCount() * iPtr->_nSpheres;

    std::vector<float> VertexColor(4*nColors, 1.0);
    
    for (size_t icol = 0; icol < nColors; ++icol)
      {
	VertexColor[4 * icol + 0] = 4.0/256.0;
	VertexColor[4 * icol + 1] = 104.0/256.0;
	VertexColor[4 * icol + 2] = 202.0/256.0;
	VertexColor[4 * icol + 3] = 1;
      }

    setGLColors(VertexColor);
  }
   
 
  {//Setup initial element data
    size_t nElements = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nElements += 3 * iPtr->_type.getFaceCount() * iPtr->_nSpheres;

    std::vector<int> ElementData(nElements, 0.0);

    nElements = 0;
    size_t nSphereVertices = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	for (size_t i = 0; i < iPtr->_nSpheres; ++i)
	  for (int j = 0; j < 3 * iPtr->_type.getFaceCount(); ++j)
	    ElementData[nElements + 3 * iPtr->_type.getFaceCount() * i + j] 
	      = i * iPtr->_type.getVertexCount() + iPtr->_type.getFaces()[j]
	      + nSphereVertices;

	nSphereVertices += iPtr->_type.getVertexCount() * iPtr->_nSpheres;
	nElements += 3 * iPtr->_type.getFaceCount() * iPtr->_nSpheres;
      }
    
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
  
  _program = cl::Program(CLState.getCommandQueue().getInfo<CL_QUEUE_CONTEXT>(), kernelSource);
  
  std::string buildOptions;
  
  cl::Device clDevice = CLState.getCommandQueue().getInfo<CL_QUEUE_DEVICE>();
  try {
    _program.build(std::vector<cl::Device>(1, clDevice), buildOptions.c_str());
  } catch(cl::Error& err) {
    
    std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLState.getDevice());
    
    std::cout << "Compilation failed for device " <<
      CLState.getDevice().getInfo<CL_DEVICE_NAME>()
	      << "\nBuild Log:" << msg;
    
    throw;
  }
  
  _renderKernel = cl::Kernel(_program, "SphereRenderKernel");
  _sortDataKernel = cl::Kernel(_program, "GenerateData");

  cl_uint paddedN = ((_N + 1023)/1024) * 1024;
  _sortDataKernelFunc = _sortDataKernel.bind(CLState.getCommandQueue(), cl::NDRange(paddedN), cl::NDRange(256));
  
  sortFunctor.build(CLState.getCommandQueue(), CLState.getContext());
  CPUsortFunctor.build(CLState.getCommandQueue(), CLState.getContext());

  for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
       iPtr != _renderDetailLevels.end(); ++iPtr)
    iPtr->setupCLBuffers(CLState);

  try {
    _sphereDataLock.unlock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially (UN)lock the sphere data.";
    }
}

void 
RTSpheres::sortTick(magnet::CL::CLGLState& CLState, const magnet::GL::viewPort& _viewPortInfo)
{
  cl_float4 campos = getclVec(_viewPortInfo._position);
  cl_float4 camdir = getclVec(_viewPortInfo._cameraDirection);
  cl_float4 camup = getclVec(_viewPortInfo._cameraUp);
  
  //Generate the sort data
  _sortDataKernelFunc(_spherePositions, _sortKeys, _sortData,
		     campos, camdir, camup,
		     (cl_float)_viewPortInfo._aspectRatio,
		     (cl_float)_viewPortInfo._zNearDist,
		     (cl_float)_viewPortInfo._fovY,
		     _N);
  
  if ((_renderDetailLevels.size() > 2) 
      || (_renderDetailLevels.front()._nSpheres != _N))
    {
      if (CLState.getCommandQueue().getInfo<CL_QUEUE_DEVICE>().getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_CPU)
	sortFunctor(_sortKeys, _sortData, _sortKeys, _sortData);
      else
	CPUsortFunctor(_sortKeys, _sortData);
    }
}

void 
RTSpheres::clTick(magnet::CL::CLGLState& CLState, const magnet::GL::viewPort&  _viewPortInfo)
{
  //First check you can get a lock on the position data!
  magnet::thread::ScopedLock lock(_sphereDataLock);

  if (!(++_frameCount % _sortFrequency)) sortTick(CLState, _viewPortInfo);
  
  clTick_no_sort_or_locking(CLState);
}

void 
RTSpheres::clTick_no_sort_or_locking(magnet::CL::CLGLState& CLState)
{
  //Aqquire GL buffer objects
  _clbuf_Positions.acquire(CLState.getCommandQueue());

  //Finally, run render kernels
  cl::KernelFunctor renderKernelFunc 
    = _renderKernel.bind(CLState.getCommandQueue(), cl::NDRange(_globalsize), cl::NDRange(_workgroupsize));

  cl_uint renderedSpheres = 0;
  cl_uint renderedVertexData = 0;
  for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
       iPtr != _renderDetailLevels.end(); ++iPtr)
    {
      cl_int vertexOffset = renderedVertexData - 3 * renderedSpheres * iPtr->_type.getVertexCount();

      renderKernelFunc(_spherePositions, (cl::Buffer)_clbuf_Positions, iPtr->_primativeVertices, 
		       iPtr->_type.getVertexCount(), renderedSpheres, renderedSpheres + iPtr->_nSpheres, 
		       vertexOffset, _sortData);

      renderedSpheres += iPtr->_nSpheres;
      renderedVertexData += 3 * iPtr->_nSpheres * iPtr->_type.getVertexCount();
    }

  //Release resources
  _clbuf_Positions.release(CLState.getCommandQueue());
}

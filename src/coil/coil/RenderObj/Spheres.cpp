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

#include "Primatives/Sphere.hpp"
#include "Spheres.clh"
#include <errno.h>

struct  SortDataType { cl_uint ID; cl_float dist;};

cl_float4 getclVec(Vector vec)
{ 
  cl_float4 clvec;
  clvec.x = vec[0];
  clvec.y = vec[1];
  clvec.z = vec[2];
  clvec.w = 0;
  return clvec;
}


RTSpheres::RTSpheres(const magnet::GL::viewPort& viewPortInfo,
		     size_t N, const std::vector<SphereDetails>& renderDetailLevels):
  _N(N),
  _renderDetailLevels(renderDetailLevels),
  _frameCount(0),
  _sortFrequency(1),
  _workgroupsize(0),
  _globalsize(0),
  _viewPortInfo(viewPortInfo)
{
  //We lock the mutex straight away until initOpenCL() has been called!
  //This stops the data being accessed before it even exists!

  try {
    _sphereDataLock.lock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially lock the sphere data.";
    }
}

void
RTSpheres::initOpenCL(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers)
{
  {
    _spherePositions = cl::Buffer(Context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, 
				  sizeof(cl_float4) *  _N);

    //We must pad the sort data out to a multiple of 1024

    size_t paddedN = ((_N + 1023)/1024) * 1024;
    
    _sortKeys = cl::Buffer(Context, CL_MEM_READ_WRITE, sizeof(float) * paddedN);
    _sortData = cl::Buffer(Context, CL_MEM_READ_WRITE, sizeof(cl_uint) * paddedN);
    
    cl_float4* Pos = (cl_float4*)CmdQ.enqueueMapBuffer(_spherePositions, true, 
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
    CmdQ.enqueueUnmapMemObject(_spherePositions, (void*)Pos);
  }

  {//Setup initial vertex positions
    size_t nVertice = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nVertice += iPtr->_type.n_vertices * iPtr->_nSpheres;

    std::vector<float> VertexPos(3 * nVertice, 0.0);
    setGLPositions(VertexPos);
    initOCLVertexBuffer(Context, hostTransfers);
  }
  
  {//Setup inital normal vectors
    size_t nNormals = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nNormals += iPtr->_type.n_vertices * iPtr->_nSpheres;

    std::vector<float> VertexNormals(3 * nNormals, 0.0);

    nNormals = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	for (size_t i = 0; i < iPtr->_nSpheres; ++i)
	  for (int j = 0; j < 3 * iPtr->_type.n_vertices; ++j)
	    VertexNormals[nNormals + 3 * iPtr->_type.n_vertices * i + j] = iPtr->_type.vertices[j];

	nNormals += 3 * iPtr->_nSpheres * iPtr->_type.n_vertices;
      }
    
    setGLNormals(VertexNormals);
  }


  {//Setup initial Colors
    size_t nColors = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      nColors += iPtr->_type.n_vertices * iPtr->_nSpheres;

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
      nElements += 3 * iPtr->_type.n_faces * iPtr->_nSpheres;

    std::vector<int> ElementData(nElements, 0.0);

    nElements = 0;
    size_t nSphereVertices = 0;
    for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	for (size_t i = 0; i < iPtr->_nSpheres; ++i)
	  for (int j = 0; j < 3 * iPtr->_type.n_faces; ++j)
	    ElementData[nElements + 3 * iPtr->_type.n_faces * i + j] 
	      = i * iPtr->_type.n_vertices + iPtr->_type.faces[j]
	      + nSphereVertices;

	nSphereVertices += iPtr->_type.n_vertices * iPtr->_nSpheres;
	nElements += 3 * iPtr->_type.n_faces * iPtr->_nSpheres;
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


  sortFunctor.build(CmdQ, Context);
  CPUsortFunctor.build(CmdQ,Context);

  for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
       iPtr != _renderDetailLevels.end(); ++iPtr)
    iPtr->setupCLBuffers(Context);

  sortTick(CmdQ, Context);

  clTick_no_sort_or_locking(CmdQ, Context);

  try {
    _sphereDataLock.unlock();
  } catch (std::exception& except)
    {
      M_throw() << "Failed to initially (UN)lock the sphere data.";
    }
}

void 
RTSpheres::sortTick(cl::CommandQueue& CmdQ, cl::Context& Context)
{

  cl_uint paddedN = ((_N + 1023)/1024) * 1024;

  cl::KernelFunctor sortDataKernelFunc 
    = _sortDataKernel.bind(CmdQ, cl::NDRange(paddedN), cl::NDRange(256));
  
  cl_float4 campos = getclVec(Vector(_viewPortInfo._cameraX, _viewPortInfo._cameraY, _viewPortInfo._cameraZ));
  cl_float4 camdir = getclVec(_viewPortInfo._cameraDirection);
  cl_float4 camup = getclVec(_viewPortInfo._cameraUp);
  
  //Generate the sort data
  sortDataKernelFunc(_spherePositions, _sortKeys, _sortData,
		     campos, camdir, camup,
		     (cl_float)_viewPortInfo._aspectRatio,
		     (cl_float)_viewPortInfo._zNearDist,
		     (cl_float)_viewPortInfo._fovY,
		     _N, paddedN);
  
  if ((_renderDetailLevels.size() > 2) || (_renderDetailLevels.front()._nSpheres != _N))
    {
      if (CmdQ.getInfo<CL_QUEUE_DEVICE>().getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_CPU)
	sortFunctor(_sortKeys, _sortData, _sortKeys, _sortData);
      else
	CPUsortFunctor(_sortKeys, _sortData);
    }
}

void 
RTSpheres::clTick(cl::CommandQueue& CmdQ, cl::Context& Context)
{
  //First check you can get a lock on the position data!
  magnet::thread::ScopedLock lock(_sphereDataLock);

  if (!(++_frameCount % _sortFrequency))
    sortTick(CmdQ, Context);
  
  clTick_no_sort_or_locking(CmdQ, Context);
}

void 
RTSpheres::clTick_no_sort_or_locking(cl::CommandQueue& CmdQ, cl::Context& Context)
{
  //Aqquire GL buffer objects
  _clbuf_Positions.acquire(CmdQ);

  //Finally, run render kernels
  cl::KernelFunctor renderKernelFunc 
    = _renderKernel.bind(CmdQ, cl::NDRange(_globalsize), cl::NDRange(_workgroupsize));

  cl_uint renderedSpheres = 0;
  cl_uint renderedVertexData = 0;
  for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
       iPtr != _renderDetailLevels.end(); ++iPtr)
    {
      cl_int vertexOffset = renderedVertexData - 3 * renderedSpheres * iPtr->_type.n_vertices;

      renderKernelFunc(_spherePositions, (cl::Buffer)_clbuf_Positions, iPtr->_primativeVertices, 
		       iPtr->_type.n_vertices, renderedSpheres, renderedSpheres + iPtr->_nSpheres, 
		       vertexOffset, _sortData);

      renderedSpheres += iPtr->_nSpheres;
      renderedVertexData += 3 * iPtr->_nSpheres * iPtr->_type.n_vertices;
    }

  //Release resources
  _clbuf_Positions.release(CmdQ);
}

cl_float4* 
RTSpheres::writePositionData(cl::CommandQueue& cmdq)
{
  _sphereDataLock.lock();
  return (cl_float4*) cmdq.enqueueMapBuffer(_spherePositions, true, CL_MAP_WRITE, 0, _N);
}

void 
RTSpheres::returnPositionData(cl::CommandQueue& cmdq, cl_float4* clBufPointer)
{
  cmdq.enqueueUnmapMemObject(_spherePositions, (void*)clBufPointer);  
  _sphereDataLock.unlock();
}

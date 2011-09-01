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

#include "Spheres.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <fstream>

#include "Spheres.clh"
#include <errno.h>
#include <coil/RenderObj/console.hpp>

namespace coil {
  struct  SortDataType { cl_uint ID; cl_float dist; };
  namespace {
    inline cl_float4 getclVec(Vector vec)
    { 
      cl_float4 clvec;
      clvec.x = vec[0];
      clvec.y = vec[1];
      clvec.z = vec[2];
      clvec.w = 0;
      return clvec;
    }
  }


  RTSpheres::RTSpheres(size_t N, std::string name):
    RTriangles(name),
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
	_renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::objects::primitives::Sphere::icosahedron, 2, stage_spheres));
	spheres_rendered += stage_spheres;
      }

    stage_spheres = std::min(size_t(1000ul), N - spheres_rendered);
    if (stage_spheres)
      {
	_renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::objects::primitives::Sphere::icosahedron, 1, stage_spheres));
	spheres_rendered += stage_spheres;
      }

    stage_spheres = std::min(size_t(10000ul), N - spheres_rendered);
    if (stage_spheres)
      {
	_renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::objects::primitives::Sphere::icosahedron, 0, stage_spheres));
	spheres_rendered += stage_spheres;
      }

    stage_spheres = std::min(size_t(200000ul), N - spheres_rendered);
    if (stage_spheres)
      {
	_renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::objects::primitives::Sphere::octahedron, 0, stage_spheres));
	spheres_rendered += stage_spheres;
      }

    stage_spheres = N - spheres_rendered;
    if (stage_spheres)
      _renderDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::objects::primitives::Sphere::tetrahedron, 0, N - spheres_rendered));
  }

  void
  RTSpheres::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RTriangles::init(systemQueue);

    magnet::GL::Context& context = magnet::GL::Context::getContext();
    //Build the sort functor now so we can grab the padding
    sortFunctor.build(context.getCLCommandQueue(), context.getCLContext());
  
    //We must pad the sort data out to a multiple of sortFunctor.padding()
    cl_uint padding = std::max(sortFunctor.padding(), size_t(1024));
    cl_uint paddedN = ((_N + padding - 1) / padding) * padding;

    {
      _spherePositions = cl::Buffer(context.getCLContext(), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, 
				    sizeof(cl_float4) *  _N);

    
      _sortKeys = cl::Buffer(context.getCLContext(), CL_MEM_READ_WRITE, sizeof(cl_float) * paddedN);
      _sortData = cl::Buffer(context.getCLContext(), CL_MEM_READ_WRITE, sizeof(cl_uint) * paddedN);
      _sphereColors = cl::Buffer(context.getCLContext(), CL_MEM_READ_ONLY, sizeof(cl_uchar4) * paddedN);

      cl_float4* Pos = (cl_float4*)context.getCLCommandQueue().enqueueMapBuffer(_spherePositions, true, 
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
      context.getCLCommandQueue().enqueueUnmapMemObject(_spherePositions, (void*)Pos);
    }

    {//Setup initial vertex positions
      size_t nVertice = 0;
      for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	   iPtr != _renderDetailLevels.end(); ++iPtr)
	nVertice += iPtr->_type.getVertexCount() * iPtr->_nSpheres;

      std::vector<float> VertexPos(3 * nVertice, 0.0);
      setGLPositions(VertexPos);
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

      std::vector<GLubyte> VertexColor(nColors * 4);
    
      for (size_t icol = 0; icol < nColors; ++icol)
	{
	  VertexColor[icol * 4 + 0] = 255;
	  VertexColor[icol * 4 + 1] = 255;
	  VertexColor[icol * 4 + 2] = 255;
	  VertexColor[icol * 4 + 3] = 255;
	}

      setGLColors(VertexColor);
    }
   
 
    {//Setup initial element data
      size_t nElements = 0;
      for (std::vector<SphereDetails>::const_iterator iPtr = _renderDetailLevels.begin();
	   iPtr != _renderDetailLevels.end(); ++iPtr)
	nElements += 3 * iPtr->_type.getFaceCount() * iPtr->_nSpheres;

      std::vector<GLuint> ElementData(nElements, 0.0);

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
  
    _program = cl::Program(context.getCLCommandQueue().getInfo<CL_QUEUE_CONTEXT>(), kernelSource);
  
    std::string buildOptions;
  
    cl::Device clDevice = context.getCLCommandQueue().getInfo<CL_QUEUE_DEVICE>();
    try {
      _program.build(std::vector<cl::Device>(1, clDevice), buildOptions.c_str());
    } catch(cl::Error& err) {
    
      std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(context.getCLDevice());
    
      std::cout << "Compilation failed for device " <<
	context.getCLDevice().getInfo<CL_DEVICE_NAME>()
		<< "\nBuild Log:" << msg;
    
      throw;
    }
  
    _renderKernel = cl::Kernel(_program, "SphereRenderKernel");
    _sortDataKernel = cl::Kernel(_program, "GenerateData");
    _colorKernel = cl::Kernel(_program, "SphereColorKernel");
    _pickingKernel = cl::Kernel(_program, "SpherePickingKernel");

    _sortDataKernelFunc = _sortDataKernel.bind(context.getCLCommandQueue(), cl::NDRange(paddedN), 
					       cl::NDRange(256));

    _renderKernelFunc = _renderKernel.bind(context.getCLCommandQueue(), cl::NDRange(_globalsize), 
					   cl::NDRange(_workgroupsize));

    _pickingKernelFunc = _pickingKernel.bind(context.getCLCommandQueue(), cl::NDRange(_globalsize), 
					     cl::NDRange(_workgroupsize));

    for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      iPtr->setupCLBuffers(context);
  }

  void 
  RTSpheres::sortTick(const magnet::GL::Camera& camera)
  {
    cl_float4 campos = getclVec(camera.getEyeLocation());
    cl_float4 camdir = getclVec(camera.getCameraDirection());
    cl_float4 camup = getclVec(camera.getCameraUp());
  
    //Generate the sort data
    _sortDataKernelFunc(_spherePositions, _sortKeys, _sortData,
			campos, camdir, camup,
			(cl_float)camera.getAspectRatio(),
			(cl_float)camera.getZNear(),
			(cl_float)camera.getFOVY(),
			_N);
  
    if ((_renderDetailLevels.size() > 2) 
	|| (_renderDetailLevels.front()._nSpheres != _N))
      sortFunctor(_sortKeys, _sortData);
  
    recolor();
  }

  void
  RTSpheres::recolor()
  {
    magnet::GL::Context& context = getContext();
    //Aqquire GL buffer objects
 
    //Run color kernels
    _colorKernelFunc = _colorKernel.bind(context.getCLCommandQueue(), 
					 cl::NDRange(_globalsize), 
					 cl::NDRange(_workgroupsize));
  
    cl::Buffer clbuff = _colBuff.acquireCLObject();

    cl_uint renderedSpheres = 0;
    cl_uint renderedVertexData = 0;
    for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	cl_int vertexOffset = renderedVertexData - renderedSpheres * iPtr->_type.getVertexCount();

	_colorKernelFunc(clbuff, _sphereColors, iPtr->_type.getVertexCount(), 
			 renderedSpheres, renderedSpheres + iPtr->_nSpheres, 
			 vertexOffset, _sortData, _N);

	renderedSpheres += iPtr->_nSpheres;
	renderedVertexData += iPtr->_nSpheres * iPtr->_type.getVertexCount();
      }

    //Release resources
    _colBuff.releaseCLObject();
  }

  void 
  RTSpheres::clTick(const magnet::GL::Camera& camera)
  {
    if (!_visible) return;

    if (!(++_frameCount % _sortFrequency)) sortTick(camera);
  
    //Aqquire GL buffer objects
    //Finally, run render kernels
    cl_uint renderedSpheres = 0;
    cl_uint renderedVertexData = 0;

    cl::Buffer clposbuf = _posBuff.acquireCLObject();

    for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	cl_int vertexOffset = renderedVertexData - 3 * renderedSpheres * iPtr->_type.getVertexCount();

	_renderKernelFunc(_spherePositions, clposbuf, iPtr->_primativeVertices, 
			  iPtr->_type.getVertexCount(), renderedSpheres, renderedSpheres + iPtr->_nSpheres, 
			  vertexOffset, _sortData);

	renderedSpheres += iPtr->_nSpheres;
	renderedVertexData += 3 * iPtr->_nSpheres * iPtr->_type.getVertexCount();
      }

    //Release resources
    _posBuff.releaseCLObject();
  }

  void 
  RTSpheres::pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, uint32_t& offset)
  {
    //Run color kernels  
    cl_uint renderedSpheres = 0;
    cl_uint renderedVertexData = 0;
  
    cl::Buffer clcolbuf = _colBuff.acquireCLObject();
    cl_uint cl_offset = offset;


    for (std::vector<SphereDetails>::iterator iPtr = _renderDetailLevels.begin();
	 iPtr != _renderDetailLevels.end(); ++iPtr)
      {
	cl_int vertexOffset = renderedVertexData - renderedSpheres * iPtr->_type.getVertexCount();
      
	_pickingKernelFunc(clcolbuf, iPtr->_type.getVertexCount(), 
			   renderedSpheres, renderedSpheres + iPtr->_nSpheres, 
			   vertexOffset, _sortData, cl_offset, _N);
      
	renderedSpheres += iPtr->_nSpheres;
	renderedVertexData += iPtr->_nSpheres * iPtr->_type.getVertexCount();
      }
  
    //Release resources
    _colBuff.releaseCLObject();
    _colBuff.getContext().getCLCommandQueue().finish();
    offset += _N;

    glRender(fbo, cam, PICKING_PASS);
  }

  void RTSpheres::finishPicking(uint32_t& offset, const uint32_t val)
  {
    recolor();

    bool picked = (val >= offset) && ((val - offset) < _N);

    if (picked)
      std::cout << "You picked a sphere! with an ID of " 
		<< (val - offset) << std::endl;

    offset += _N;
  }
}

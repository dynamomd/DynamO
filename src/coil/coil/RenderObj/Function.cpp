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

#include "Function.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

namespace coil {
  RFunction::RFunction(size_t N, Vector origin, Vector axis1,
		       Vector axis2, Vector axis3, cl_float functionOriginX,
		       cl_float functionOriginY, cl_float functionRangeX, 
		       cl_float functionRangeY, bool drawAxis, bool staticShape,
		       std::string name,
		       std::string function,
		       std::string normalCalc,
		       std::string colorCalc):
    RTriangles(name),
    _origin(origin),
    _axis1(axis1),
    _axis2(axis2),
    _axis3(axis3),
    _drawAxis(drawAxis),
    _staticShape(staticShape),
    _function(function),
    _normalCalc(normalCalc),
    _colorCalc(colorCalc)    
  {
    //Copy to the cl types
    for (size_t i(0); i < 3; ++i)
      {
	_cl_origin.s[i] = origin[i];
	_cl_axis1.s[i] = axis1[i];
	_cl_axis2.s[i] = axis2[i];
	_cl_axis3.s[i] = axis3[i];
      }

    _cl_origin.s[3] = 0;
    _cl_axis1.s[3]  = 0;
    _cl_axis2.s[3]  = 0;
    _cl_axis3.s[3]  = 0;

    _functionOrigin.s[0] = functionOriginX;
    _functionOrigin.s[1] = functionOriginY;
    _functionRange.s[0] = functionRangeX;
    _functionRange.s[1] = functionRangeY;

    //Set N to a multiple of 16 so the workgroup size of 256 always fits
    _N = (N + 15) / 16;
    _N *= 16;
  }

  void 
  RFunction::initOpenGL()
  {}

  void 
  RFunction::initOpenCL()
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
      std::vector<GLubyte> VertexColor(_N * _N * 4);
      for (size_t i = 0; i < _N; i++)
	{       
	  for (size_t j = 0; j < _N; j++)
	    {
	      VertexColor[(i + _N * j) * 4 + 0] = 255;
	      VertexColor[(i + _N * j) * 4 + 1] = 255;
	      VertexColor[(i + _N * j) * 4 + 2] = 255;
	      VertexColor[(i + _N * j) * 4 + 3] = 255;
	    }       
	}
      setGLColors(VertexColor);
    }
   
    {//Setup initial element data
      std::vector<GLuint> ElementData(3*2*(_N-1)*(_N-1), 0.0);
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
  
    _kernelsrc = genKernelSrc();

    cl::Program::Sources kernelSource;

    kernelSource.push_back(std::pair<const char*, ::size_t>(_kernelsrc.c_str(), _kernelsrc.size()));
  
    _program = cl::Program(magnet::GL::Context::getContext().getCLContext(), kernelSource);

    try
      { _program.build(std::vector<cl::Device>(1, magnet::GL::Context::getContext().getCLDevice())); }
    catch(cl::Error& err) 
      {    
	std::string msg 
	  = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(magnet::GL::Context::getContext().getCLDevice());
      
	std::cout << "Compilation failed for device " <<
	  magnet::GL::Context::getContext().getCLDevice().getInfo<CL_DEVICE_NAME>()
		  << "\nBuild Log:" << msg;
	throw;
      }
  
    _kernel = cl::Kernel(_program, "FunctionRenderKernel");
    _pickKernel = cl::Kernel(_program, "FunctionPickKernel");
    const size_t workgroupSize = 256;
  
    magnet::GL::Context& context = magnet::GL::Context::getContext();
    //N is a multiple of 16, so a workgroup size of 256 is always good
    _kernelFunc = _kernel.bind(context.getCLCommandQueue(), cl::NDRange(_N * _N), 
			       cl::NDRange(workgroupSize));

    _pickFunc = _pickKernel.bind(context.getCLCommandQueue(), cl::NDRange(_N * _N), 
				 cl::NDRange(workgroupSize));
  
    clock_gettime(CLOCK_MONOTONIC, &startTime);
    glFinish();
  
    //Now do the first clTick if the shape is static!
    if (_staticShape)
      {
	_staticShape = false;
	clTick();
	_staticShape = true;      
      }
  }


  void 
  RFunction::clTick()
  {
    if (_staticShape || !_visible) return;

    timespec currTime;
    clock_gettime(CLOCK_MONOTONIC, &currTime);
  
    tempo = float(currTime.tv_sec) - float(startTime.tv_sec)
      + 1e-9 * (float(currTime.tv_nsec) - float(startTime.tv_nsec));

    //Run Kernel
    _kernelFunc((cl::Buffer)_posBuff.acquireCLObject(), 
		(cl::Buffer)_colBuff.acquireCLObject(), 
		(cl::Buffer)_normBuff.acquireCLObject(), 
		tempo,
		_functionOrigin,
		_functionRange,
		_cl_axis1,
		_cl_axis2,
		_cl_axis3,
		_cl_origin,
		_N, _A);
  
    //Release resources
    _posBuff.releaseCLObject();
    _colBuff.releaseCLObject();
    _normBuff.releaseCLObject();
  }

  void 
  RFunction::glRender()
  {
    RTriangles::glRender();

    //Draw the axis of the plotted function
    if (_drawAxis && _visible)
      {
	coil::glprimatives::drawArrow(_origin, _origin + _axis1);
	coil::glprimatives::drawArrow(_origin, _origin + _axis2);
	coil::glprimatives::drawArrow(_origin, _origin + _axis3);
      }
  }

  std::string
  RFunction::genKernelSrc()
  {

#define STRINGIFY(A) #A

    return std::string(STRINGIFY(

				 __kernel void
				 FunctionRenderKernel(__global float * positions,
						      __global uchar4 * colors,
						      __global float * normals,
						      float t,
						      float2 functionOrigin,
						      float2 functionRange,
						      float4 axis1,
						      float4 axis2,
						      float4 axis3,
						      float4 origin,
						      uint N, float A)
				 {
				   positions += 3 * get_global_id(0);
				   normals += 3 * get_global_id(0);
				   colors += get_global_id(0);

				   float2 normPos = (float2)(get_global_id(0) % N, get_global_id(0) / N);
				   normPos /= N;

				   float2 pos = normPos * functionRange + functionOrigin;

				   float f; 
				   )) + _function + std::string(STRINGIFY(
  
									  float4 vertexPosition = normPos.x * axis1 + normPos.y * axis2 + f * axis3 + origin;

									  positions[0] = vertexPosition.x;
									  positions[1] = vertexPosition.y;
									  positions[2] = vertexPosition.z;

									  float4 normal;
									  )) + _normalCalc + std::string(STRINGIFY(
														   normal *= (float4)(functionRange * length(axis3) , 1.0f / length(axis3), 0);

														   float4 rotatedNormal 
														   = normalize(normal.x * axis1 +
															       normal.y * axis2 +
															       normal.z * axis3
															       );

														   normals[0] = rotatedNormal.x;
														   normals[1] = rotatedNormal.y;
														   normals[2] = rotatedNormal.z;

														   )) + _colorCalc + std::string(STRINGIFY(
																			   }

																		 __kernel void
																		 FunctionPickKernel(__global uint * colors, uint offset)
																		 {
																		   colors[get_global_id(0)] = get_global_id(0) + offset;
																		 }
																		 ));
  }


  void 
  RFunction::initPicking(cl_uint& offset)
  {
    //Run Kernel
    _pickFunc((cl::Buffer)_colBuff.acquireCLObject(), offset);
    //Release resources
    _colBuff.releaseCLObject();

    offset += _N * _N;
  }

  void 
  RFunction::pickingRender()
  {
    RTriangles::glRender();
  }

  void 
  RFunction::finishPicking(cl_uint& offset, const cl_uint val)
  {
    if (val - offset < _N * _N)
      std::cout << "You picked a function point with coords of " 
		<< (val - offset) % _N
		<< ","
		<< (val - offset) / _N 
		<< std::endl;

    //Run Kernel
    _kernelFunc((cl::Buffer)_posBuff.acquireCLObject(),
		(cl::Buffer)_colBuff.acquireCLObject(),
		(cl::Buffer)_normBuff.acquireCLObject(),
		tempo,
		_functionOrigin,
		_functionRange,
		_cl_axis1,
		_cl_axis2,
		_cl_axis3,
		_cl_origin,
		_N);
  
    //Release resources
    _posBuff.releaseCLObject();
    _colBuff.releaseCLObject();
    _normBuff.releaseCLObject();

    offset += _N * _N;
  }
}

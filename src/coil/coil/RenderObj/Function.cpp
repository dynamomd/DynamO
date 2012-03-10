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
#include <magnet/string/formatcode.hpp>
#include <magnet/string/line_number.hpp>


extern const guint8 Function_Icon[];
extern const size_t Function_Icon_size;

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RFunction::getIcon()
  {
    return Gdk::Pixbuf::create_from_inline(Function_Icon_size, Function_Icon);
  }


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
    _program(function, normalCalc, colorCalc)
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
  RFunction::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RTriangles::init(systemQueue);

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
    magnet::GL::Context::ContextPtr context = magnet::GL::Context::getContext();

    _program.build(context->getCLCommandQueue(), 
		   context->getCLContext());
    _kernel = _program["FunctionRenderKernel"];
    _pickKernel = _program["FunctionPickKernel"];
    const size_t workgroupSize = 256;
  
    //N is a multiple of 16, so a workgroup size of 256 is always good
    _kernelFunc = _kernel.bind(context->getCLCommandQueue(), cl::NDRange(_N * _N), 
			       cl::NDRange(workgroupSize));

    _pickFunc = _pickKernel.bind(context->getCLCommandQueue(), cl::NDRange(_N * _N), 
				 cl::NDRange(workgroupSize));
  
    clock_gettime(CLOCK_MONOTONIC, &startTime);
    glFinish();
  
    //Now do the first clTick if the shape is static!
    if (_staticShape)
      {
	_staticShape = false;
	//Ignore the state of the visible flag
	bool isVisible = _visible;
	_visible = true;
	clTick();
	_visible = isVisible;
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
  RFunction::glRender(const magnet::GL::Camera& cam, RenderMode mode)
  {
    RTriangles::glRender(cam, mode);

    //Draw the axis of the plotted function
    if (_drawAxis && _visible)
      {
	coil::glprimatives::drawArrow(_origin, _origin + _axis1);
	coil::glprimatives::drawArrow(_origin, _origin + _axis2);
	coil::glprimatives::drawArrow(_origin, _origin + _axis3);
      }
  }

//  void 
//  RFunction::pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, uint32_t& offset)
//  {
//    //Run Kernel
//    glFinish();
//    cl_uint cl_offset = offset;
//    _pickFunc((cl::Buffer)_colBuff.acquireCLObject(), cl_offset);
//    //Release resources
//    _colBuff.releaseCLObject();
//    offset += _N * _N;
//    _colBuff.getContext().getCLCommandQueue().finish();
//
//    RTriangles::glRender(fbo, cam, PICKING_PASS);
//  }
//
//  bool
//  RFunction::finishPicking(uint32_t& offset, const uint32_t val)
//  {
//    //Run Kernel
//    _kernelFunc((cl::Buffer)_posBuff.acquireCLObject(),
//		(cl::Buffer)_colBuff.acquireCLObject(),
//		(cl::Buffer)_normBuff.acquireCLObject(),
//		tempo,
//		_functionOrigin,
//		_functionRange,
//		_cl_axis1,
//		_cl_axis2,
//		_cl_axis3,
//		_cl_origin,
//		_N);
//  
//    //Release resources
//    _posBuff.releaseCLObject();
//    _colBuff.releaseCLObject();
//    _normBuff.releaseCLObject();
//
//    bool picked = (val >= offset) && ((val - offset) < (_N * _N));
//
//    if (picked)
//      std::cout << "You picked a function point (id=" << (val - offset)
//		<< ") with coords of " << (val - offset) % _N
//		<< "," << (val - offset) / _N 
//		<< std::endl;
//    
//    offset += _N * _N;
//    
//    return picked;
//  }
}

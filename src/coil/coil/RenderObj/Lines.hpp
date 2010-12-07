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
#pragma once
#include "RenderObj.hpp"
#include <vector>

#include <magnet/CL/GLBuffer.hpp>

class RLines : public RenderObj
{
public:
  RLines(size_t N);
  ~RLines();

  virtual void glRender();

  virtual void initOpenGL();

  virtual void initOpenCL(magnet::CL::CLGLState&) {}
  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&) {}

  void setGLColors(std::vector<cl_uchar4>& VertexColor);
  void setGLPositions(std::vector<float>& VertexPos);
  void setGLElements(std::vector<int>& Elements);

  void initOCLVertexBuffer(cl::Context Context);
  void initOCLColorBuffer(cl::Context Context);
  void initOCLElementBuffer(cl::Context Context);

protected:
  size_t _N;
  GLuint _colBuff;
  size_t _colBuffSize;
  cl::GLBuffer _clbuf_Colors;

  GLuint _posBuff;
  size_t _posBuffSize;
  cl::GLBuffer _clbuf_Positions;
  
  GLuint _elementBuff;
  size_t _elementBuffSize;
  cl::GLBuffer _clbuf_Elements;
};

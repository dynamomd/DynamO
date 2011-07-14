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
#pragma once
#include "RenderObj.hpp"
#include <vector>

#include <magnet/CL/GLBuffer.hpp>

class RLines : public RenderObj
{
public:
  RLines(size_t N, std::string name);
  ~RLines();

  virtual void glRender();
  virtual void initOpenGL();
  virtual void initOpenCL();

  void setGLColors(std::vector<cl_uchar4>& VertexColor);
  void setGLPositions(std::vector<GLfloat>& VertexPos);
  void setGLElements(std::vector<GLuint>& Elements);

  void initOCLVertexBuffer(cl::Context Context);
  void initOCLColorBuffer(cl::Context Context);
  void initOCLElementBuffer(cl::Context Context);

  magnet::GL::Buffer<GLfloat>& getVertexGLData() { return _posBuff; }
  magnet::GL::Buffer<cl_uchar4>& getColorGLData() { return _colBuff; }

  virtual void releaseCLGLResources();

protected:
  size_t _N;
  magnet::GL::Buffer<cl_uchar4> _colBuff;
  cl::GLBuffer<cl_uchar4> _clbuf_Colors;
  
  magnet::GL::Buffer<GLfloat> _posBuff;
  cl::GLBuffer<GLfloat> _clbuf_Positions;
  
  magnet::GL::Buffer<GLuint> _elementBuff;
  cl::GLBuffer<GLuint> _clbuf_Elements;
};

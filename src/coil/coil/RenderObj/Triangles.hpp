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
#include <gtkmm.h>
#include "RenderObj.hpp"
#include <magnet/GL/buffer.hpp>
#include <vector>
#include <memory>
#include <magnet/CL/GLBuffer.hpp>

class RTriangles : public RenderObj
{
public:
  RTriangles(std::string name);
  ~RTriangles();

  virtual void glRender();

  void setGLColors(std::vector<cl_uchar4>& VertexColor);
  void setGLPositions(std::vector<float>& VertexPos);
  void setGLNormals(std::vector<float>& VertexNormals);
  void setGLElements(std::vector<GLuint>& Elements);

  void initOCLVertexBuffer(cl::Context Context);
  void initOCLColorBuffer(cl::Context Context);
  void initOCLNormBuffer(cl::Context Context);
  void initOCLElementBuffer(cl::Context Context);

  virtual void releaseCLGLResources();

  virtual void initGTK();

  virtual void showControls(Gtk::ScrolledWindow* win);

  virtual void setRenderMode(RenderModeType rm);

  virtual void initPicking(cl_uint& offset);
  virtual void pickingRender();
  virtual void finishPicking(cl_uint& offset, const cl_uint val);

protected:
  void guiUpdate();
  
  std::auto_ptr<Gtk::VBox>        _gtkOptList;
  std::auto_ptr<Gtk::RadioButton> _gtkLineRender;
  std::auto_ptr<Gtk::RadioButton> _gtkPointRender;
  std::auto_ptr<Gtk::RadioButton> _gtkTriangleRender;

  magnet::GL::Buffer<cl_uchar4> _colBuff;
  cl::GLBuffer<cl_uchar4> _clbuf_Colors;

  magnet::GL::Buffer<GLfloat> _posBuff;
  cl::GLBuffer<GLfloat> _clbuf_Positions;
  
  magnet::GL::Buffer<GLfloat> _normBuff;
  cl::GLBuffer<GLfloat> _clbuf_Normals;

  magnet::GL::Buffer<GLuint> _elementBuff;
  cl::GLBuffer<GLuint> _clbuf_Elements;

  magnet::GL::Buffer<GLuint> _specialElementBuff;

  bool _pickingRenderMode;
  magnet::GL::Buffer<cl_uchar4> _pickingColorBuff;

};

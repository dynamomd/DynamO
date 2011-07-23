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

namespace coil {
  class RTriangles : public RenderObj
  {
  public:
    RTriangles(std::string name);
    ~RTriangles();

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam);
    virtual void clTick(const magnet::GL::Camera& cam) {}

    void setGLColors(std::vector<GLubyte>& VertexColor);
    void setGLPositions(std::vector<float>& VertexPos);
    void setGLNormals(std::vector<float>& VertexNormals);
    void setGLElements(std::vector<GLuint>& Elements);

    virtual void releaseCLGLResources();

    virtual void initGTK();

    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual void setRenderMode(RenderModeType rm);

    virtual void initPicking(cl_uint& offset);
    virtual void pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam);
    virtual void finishPicking(cl_uint& offset, const cl_uint val);

    magnet::GL::Context& getContext() { return _posBuff.getContext(); }

  protected:
    void guiUpdate();
  
    std::auto_ptr<Gtk::VBox>        _gtkOptList;
    std::auto_ptr<Gtk::RadioButton> _gtkLineRender;
    std::auto_ptr<Gtk::RadioButton> _gtkPointRender;
    std::auto_ptr<Gtk::RadioButton> _gtkTriangleRender;

    magnet::GL::Buffer<GLubyte> _colBuff;
    magnet::GL::Buffer<GLfloat> _posBuff;
    magnet::GL::Buffer<GLfloat> _normBuff;
    magnet::GL::Buffer<GLuint> _elementBuff;
    magnet::GL::Buffer<GLuint> _specialElementBuff;

    bool _pickingRenderMode;
    magnet::GL::Buffer<GLubyte> _pickingColorBuff;

  };
}


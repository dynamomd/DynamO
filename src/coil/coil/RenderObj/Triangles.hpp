/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <magnet/GL/shader/render.hpp>
#include <vector>
#include <memory>

namespace coil {
  class RTriangles : public RenderObj
  {
  public:
    RTriangles(std::string name);
    ~RTriangles();

    virtual void glRender(const magnet::GL::Camera& cam, RenderMode mode, const uint32_t offset);

    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    void setGLColors(const std::vector<GLubyte>& VertexColor);
    void setGLPositions(const std::vector<float>& VertexPos);
    void setGLNormals(const std::vector<float>& VertexNormals);
    void setGLElements(const std::vector<GLuint>& Elements);

    virtual void deinit();

    void setComponents(size_t comps) { _triangleComponents = comps; }

    virtual void showControls(Gtk::ScrolledWindow* win);

    const magnet::GL::Context::ContextPtr& getContext() const { return _posBuff.getContext(); }

  protected:
    enum RenderModeType { POINTS, LINES, TRIANGLES };
    RenderModeType _RenderMode;

    void setRenderMode(RenderModeType rm);

    void initGTK();
    void guiUpdate();
  
    std::unique_ptr<Gtk::VBox>        _gtkOptList;
    std::unique_ptr<Gtk::RadioButton> _gtkLineRender;
    std::unique_ptr<Gtk::RadioButton> _gtkPointRender;
    std::unique_ptr<Gtk::RadioButton> _gtkTriangleRender;

    magnet::GL::Buffer<GLubyte> _colBuff;
    magnet::GL::Buffer<GLfloat> _posBuff;
    magnet::GL::Buffer<GLfloat> _normBuff;
    magnet::GL::Buffer<GLuint> _elementBuff;
    magnet::GL::Buffer<GLuint> _specialElementBuff;

    magnet::GL::shader::RenderShader _renderShader;
    magnet::GL::shader::RenderVSMShader _renderVSMShader;
    size_t _triangleComponents;
  };
}


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
#include <magnet/GL/buffer.hpp>
#include <magnet/GL/shader/render.hpp>
#include <coil/RenderObj/RenderObj.hpp>
#include <array>
#include <memory>
#include <sstream>
#include <list>
#include <iostream>

namespace coil {
  class Console: public RenderObj
  {
  public:
    struct end {};
    
    inline Console(magnet::GL::Context::ContextPtr context, std::array<GLfloat, 3> color):
      RenderObj(context, "Console")
    {}
    
    void interfaceRender(const magnet::GL::Camera& cam, magnet::GL::objects::CairoSurface& cairo);

    void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    void showControls(Gtk::ScrolledWindow* win);
    void deinit();
    void glRender(const magnet::GL::Camera& cam, RenderMode mode, const uint32_t offset = 0);
    
  private:
    void initGTK();
    void guiUpdate();

    int _glutLastTime;

    magnet::GL::shader::RenderShader _renderShader;
    magnet::GL::Buffer<GLfloat> _gridVertices;
    std::unique_ptr<Gtk::VBox> _optList;
    std::unique_ptr<Gtk::CheckButton> _showGrid;
    std::unique_ptr<Gtk::CheckButton> _showAxis;
  };
}

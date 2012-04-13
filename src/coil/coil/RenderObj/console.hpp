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
#include <coil/RenderObj/RenderObj.hpp>
#include <magnet/GL/objects/axis.hpp>
#include <magnet/GL/objects/grid.hpp>
#include <magnet/GL/objects/cairo.hpp>
#include <tr1/array>
#include <memory>
#include <sstream>
#include <list>
#include <iostream>

namespace coil {
  class Console: public RenderObj
  {

    class AxisText : public magnet::GL::objects::CairoSurface
    {
    public:
      inline void init(size_t width = 100, size_t height = 100, size_t alpha_testing = 0) 
      { CairoSurface::init(100,100,0); }

      inline void glRender(const magnet::GL::GLMatrix& viewProjection)
      {
	{
	  std::tr1::array<GLfloat, 4> vec = {{0.55,-0.4,-0.4,1.0}};
	  vec = viewProjection * vec;
	  Xx = 0.5 + 0.5 * vec[0] / vec[3];
	  Xy = 0.5 - 0.5 * vec[1] / vec[3];
	}
	{
	  std::tr1::array<GLfloat, 4> vec = {{-0.4,0.55,-0.4,1.0}};
	  vec = viewProjection * vec;
	  Yx = 0.5 + 0.5 * vec[0] / vec[3];
	  Yy = 0.5 - 0.5 * vec[1] / vec[3];
	}
	{
	  std::tr1::array<GLfloat, 4> vec = {{-0.4,-0.4,0.55,1.0}};
	  vec = viewProjection * vec;
	  Zx = 0.5 + 0.5 * vec[0] / vec[3];
	  Zy = 0.5 - 0.5 * vec[1] / vec[3];
	}
	CairoSurface::redraw();
	CairoSurface::glRender();
      }

    protected:
      virtual void drawCommands() 
      {
	_cairoContext->scale(_width,_height);
	_cairoContext->set_font_size(0.2);

	_cairoContext->set_source_rgba(1.0, 0.1, 0.1, 1);
	_cairoContext->move_to(Xx,Xy);
	_cairoContext->show_text("X");
	_cairoContext->set_source_rgba(0.1, 1.0, 0.1, 1);
	_cairoContext->move_to(Yx,Yy);
	_cairoContext->show_text("Y");
	_cairoContext->set_source_rgba(0.1, 0.1, 1.0, 1);
	_cairoContext->move_to(Zx,Zy);
	_cairoContext->show_text("Z");
      }

      GLfloat Xx,Xy,Yx,Yy,Zx,Zy;
    };

  public:
    struct end {};

    inline Console(std::tr1::array<GLfloat, 3> color):
      RenderObj("Console")
    {}
    
    void interfaceRender(const magnet::GL::Camera& cam);

    void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    void showControls(Gtk::ScrolledWindow* win);
    void deinit() { _axis.deinit(); _grid.deinit(); }
    void glRender(const magnet::GL::Camera& cam, RenderMode mode);
    
  private:
    void initGTK();
    void guiUpdate();

    int _glutLastTime;

    magnet::GL::objects::Axis _axis;
    magnet::GL::objects::Grid _grid;
    AxisText _cairoOverlay;

    std::auto_ptr<Gtk::VBox> _optList; 
    std::auto_ptr<Gtk::CheckButton> _showGrid;
    std::auto_ptr<Gtk::CheckButton> _showAxis;
  };
}

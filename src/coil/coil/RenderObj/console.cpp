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

#include <coil/RenderObj/console.hpp>
#include <magnet/exception.hpp>
#include <magnet/clamp.hpp>
#include <magnet/GL/objects/cairo.hpp>

namespace coil {
  void 
  Console::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue);
    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);
    initGTK();
    _renderShader.build();
    _initialised = true;
  }

  namespace{
    void local_move_to(magnet::GL::objects::CairoSurface& cairo, 
		       magnet::GL::GLMatrix projViewMatrix, 
		       float x, float y, float z)
      {
	magnet::math::NVector<GLfloat, 4> vec = {x,y,z,1.0};
	vec = projViewMatrix * vec;
	cairo.getContext()->move_to(0.5 + 0.5 * vec[0] / vec[3], 0.5 - 0.5 * vec[1] / vec[3]);
      }

    void local_line_to(magnet::GL::objects::CairoSurface& cairo, 
		       magnet::GL::GLMatrix projViewMatrix, 
		       float x, float y, float z)
      {
	magnet::math::NVector<GLfloat, 4> vec = {x,y,z,1.0};
	vec = projViewMatrix * vec;
	cairo.getContext()->line_to(0.5 + 0.5 * vec[0] / vec[3],
				   0.5 - 0.5 * vec[1] / vec[3]);
      }
  }

  void
  Console::deinit()
  {
    _renderShader.deinit();
  }

  void 
  Console::interfaceRender(const magnet::GL::Camera& camera, magnet::GL::objects::CairoSurface& cairo)
  {
    //Only draw if the console has something in it or if it's visible
    if (!_visible) return;

    if (_showAxis->get_active())
      {
	const GLdouble nearPlane = 0.1,
	  axisScale = 0.09;

	cairo.getContext()->save();
	  
	magnet::GL::GLMatrix projViewMatrix
	  = magnet::GL::perspective(45, 1, nearPlane, 1000)
	  * magnet::GL::translate(0, 0, -(nearPlane + axisScale))
	  * camera.getViewRotationMatrix()
	  * magnet::GL::scale(axisScale, axisScale, axisScale);

	//Scale to a 100x100 box
	cairo.getContext()->translate(0, cairo.getHeight()-100);
	cairo.getContext()->scale(100, 100);

	cairo.getContext()->set_line_width(0.015);
	cairo.getContext()->set_font_size(0.2);

	local_move_to(cairo, projViewMatrix, -0.5, -0.5, -0.5);
	local_line_to(cairo, projViewMatrix, +0.5, -0.5, -0.5);
	local_line_to(cairo, projViewMatrix, +0.25, -0.25, -0.5);
	local_move_to(cairo, projViewMatrix, +0.5, -0.5, -0.5);
	local_line_to(cairo, projViewMatrix, +0.25, -0.5, -0.25);
	cairo.getContext()->set_source_rgba(1, 0.3, 0.3, 1);
	cairo.getContext()->stroke();

	local_move_to(cairo, projViewMatrix, -0.5, -0.5, -0.5);
	local_line_to(cairo, projViewMatrix, -0.5, +0.5, -0.5);
	local_line_to(cairo, projViewMatrix, -0.25, +0.25, -0.5);
	local_move_to(cairo, projViewMatrix, -0.5, +0.5, -0.5);
	local_line_to(cairo, projViewMatrix, -0.5, +0.25, -0.25);
	cairo.getContext()->set_source_rgba(0.3, 1, 0.3, 1);
	cairo.getContext()->stroke();

	local_move_to(cairo, projViewMatrix, -0.5, -0.5, -0.5);
	local_line_to(cairo, projViewMatrix, -0.5, -0.5, +0.5);
	local_line_to(cairo, projViewMatrix, -0.25, -0.5, +0.25);
	local_move_to(cairo, projViewMatrix, -0.5, -0.5, +0.5);
	local_line_to(cairo, projViewMatrix, -0.5, -0.25, +0.25);
	cairo.getContext()->set_source_rgba(0.3, 0.3, 1, 1);
	cairo.getContext()->stroke();

	cairo.getContext()->set_source_rgba(0.3, 0.3, 0.3, 1);
	local_move_to(cairo, projViewMatrix, +0.5, -0.35, -0.35);
	cairo.getContext()->show_text("X");
	local_move_to(cairo, projViewMatrix, -0.35, +0.5, -0.35);
	cairo.getContext()->show_text("Y");
	local_move_to(cairo, projViewMatrix, -0.35, -0.35, +0.5);
	cairo.getContext()->show_text("Z");

	cairo.getContext()->restore();
      }
  }

  void Console::glRender(const magnet::GL::Camera& camera, RenderMode mode)
  {
  }


  void
  Console::initGTK()
  {
    _optList.reset(new Gtk::VBox);//The Vbox of options   
    {
      _showAxis.reset(new Gtk::CheckButton("Show axis"));
      _showAxis->set_active(true);
      _optList->pack_start(*_showAxis, false, false); 
      _showAxis->show();
    }
    
    _optList->show();
    guiUpdate();
  }

  void
  Console::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _optList->unparent();
    win->add(*_optList);
    win->show();
  }

  void 
  Console::guiUpdate()
  {}

}

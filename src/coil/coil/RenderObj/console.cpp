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

#include <coil/RenderObj/console.hpp>
#include <magnet/exception.hpp>
#include <magnet/clamp.hpp>
#include <magnet/GL/objects/cairo.hpp>
#include <coil/glprimatives/arrow.hpp>

extern const unsigned char _binary_coilfont_ttf_start[];
extern const unsigned char _binary_coilfont_ttf_end[];

namespace coil {
  void 
  Console::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue);
    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);

    _grid.init(10,10);
    initGTK();
  }

  namespace{
    void local_move_to(magnet::GL::objects::CairoSurface& cairo, 
		       magnet::GL::GLMatrix projViewMatrix, 
		       float x, float y, float z)
      {
	std::tr1::array<GLfloat, 4> vec = {{x,y,z,1.0}};
	vec = projViewMatrix * vec;
	cairo.getContext()->move_to(0.5 + 0.5 * vec[0] / vec[3],
				   0.5 - 0.5 * vec[1] / vec[3]);
      }

    void local_line_to(magnet::GL::objects::CairoSurface& cairo, 
		       magnet::GL::GLMatrix projViewMatrix, 
		       float x, float y, float z)
      {
	std::tr1::array<GLfloat, 4> vec = {{x,y,z,1.0}};
	vec = projViewMatrix * vec;
	cairo.getContext()->line_to(0.5 + 0.5 * vec[0] / vec[3],
				   0.5 - 0.5 * vec[1] / vec[3]);
      }
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
	  = magnet::GL::GLMatrix::perspective(45, 1, nearPlane, 1000)
	  * magnet::GL::GLMatrix::translate(0, 0, -(nearPlane + axisScale))
	  * camera.getViewRotationMatrix()
	  * magnet::GL::GLMatrix::scale(axisScale, axisScale, axisScale);

	//Scale to a 100x100 box
	cairo.getContext()->translate(0, camera.getHeight()-100);
	cairo.getContext()->scale(100, 100);

	cairo.getContext()->rectangle(0,0,1,1);
	cairo.getContext()->set_source_rgba(0.5, 0.5, 0.5, 0.5);
	cairo.getContext()->fill();

	cairo.getContext()->set_line_width(0.02);
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

	cairo.getContext()->set_source_rgba(1, 1, 1, 1);
	local_move_to(cairo, projViewMatrix, +0.5, -0.5, -0.5);
	cairo.getContext()->show_text("X");
	local_move_to(cairo, projViewMatrix, -0.5, +0.5, -0.5);
	cairo.getContext()->show_text("Y");
	local_move_to(cairo, projViewMatrix, -0.5, -0.5, +0.5);
	cairo.getContext()->show_text("Z");

	cairo.getContext()->restore();
      }    
  }

  void Console::glRender(const magnet::GL::Camera& camera, RenderMode mode)
  {
    if (_showGrid->get_active())
      {
	using namespace magnet::GL;
	const Context::ContextPtr& context = magnet::GL::Context::getContext();

	GLMatrix old_model_view
	  = context->getAttachedShader()["ViewMatrix"].as<GLMatrix>();

	context->getAttachedShader()["ViewMatrix"] = camera.getViewPlaneMatrix();

	context->color(1,1,1,1);
	//Back face
	context->setAttribute(Context::instanceOriginAttrIndex, 0, 0, -camera.getScreenPlaneWidth(), 0);
	context->setAttribute(Context::instanceScaleAttrIndex,
			      camera.getScreenPlaneWidth(),
			      camera.getScreenPlaneHeight(), 1);
	_grid.glRender();

	//Sides
	context->setAttribute(Context::instanceOriginAttrIndex, 
			      0.5 * camera.getScreenPlaneWidth(), 0, 
			      -0.5 * camera.getScreenPlaneWidth(), 0);
	context->rotation(M_PI / 2, Vector(0, 1, 0));
  	_grid.glRender(); //Right side
	context->setAttribute(Context::instanceOriginAttrIndex, 
			      -0.5 * camera.getScreenPlaneWidth(), 0, 
			      -0.5 * camera.getScreenPlaneWidth(), 0);
  	_grid.glRender(); //Left side

	//Top and bottom
	context->rotation(M_PI / 2, Vector(1, 0, 0));
	context->setAttribute(Context::instanceScaleAttrIndex,
			      camera.getScreenPlaneWidth(),
			      camera.getScreenPlaneWidth(), 1);
	context->setAttribute(Context::instanceOriginAttrIndex, 0,
			      -0.5 * camera.getScreenPlaneHeight(), 
			      -0.5 * camera.getScreenPlaneWidth(), 0);
	_grid.glRender();//bottom
	context->setAttribute(Context::instanceOriginAttrIndex, 0,
			      0.5 * camera.getScreenPlaneHeight(), 
			      -0.5 * camera.getScreenPlaneWidth(), 0);
	_grid.glRender();//top
	context->getAttachedShader()["ViewMatrix"] = old_model_view;
      }
  }


  void
  Console::initGTK()
  {
    _optList.reset(new Gtk::VBox);//The Vbox of options   
    {
      _showGrid.reset(new Gtk::CheckButton("Show viewing grid"));
      _showGrid->set_active(false);
      _optList->pack_start(*_showGrid,false,false); 
      _showGrid->show();
    }

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

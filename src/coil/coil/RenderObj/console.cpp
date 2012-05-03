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

    _axis.init();
    _grid.init(10,10);
    initGTK();
  }

  void 
  Console::interfaceRender(const magnet::GL::Camera& camera, magnet::GL::objects::CairoSurface& cairo)
  {
    //Only draw if the console has something in it or if it's visible
    if (!_visible) return;

    using namespace magnet::GL;
    const Context::ContextPtr& context = _axis.getContext();
    //Draw the console in orthograpic projection
    context->cleanupAttributeArrays();

    if (_showAxis->get_active())
      {
	//Lets try to render an Axis with the rotation of the

//	    class AxisText : public magnet::GL::objects::CairoSurface
//    {
//    public:
//      inline void init(size_t width = 100, size_t height = 100, size_t alpha_testing = 0) 
//      { CairoSurface::init(100,100,0); }
//
//      inline void glRender(const magnet::GL::GLMatrix& viewProjection)
//      {
//	{
//	  std::tr1::array<GLfloat, 4> vec = {{0.55,-0.4,-0.4,1.0}};
//	  vec = viewProjection * vec;
//	  Xx = 0.5 + 0.5 * vec[0] / vec[3];
//	  Xy = 0.5 - 0.5 * vec[1] / vec[3];
//	}
//	{
//	  std::tr1::array<GLfloat, 4> vec = {{-0.4,0.55,-0.4,1.0}};
//	  vec = viewProjection * vec;
//	  Yx = 0.5 + 0.5 * vec[0] / vec[3];
//	  Yy = 0.5 - 0.5 * vec[1] / vec[3];
//	}
//	{
//	  std::tr1::array<GLfloat, 4> vec = {{-0.4,-0.4,0.55,1.0}};
//	  vec = viewProjection * vec;
//	  Zx = 0.5 + 0.5 * vec[0] / vec[3];
//	  Zy = 0.5 - 0.5 * vec[1] / vec[3];
//	}
//
//	CairoSurface::clear();
//	_cairoContext->save();
//	_cairoContext->scale(_width,_height);
//	_cairoContext->set_font_size(0.2);
//
//	_cairoContext->set_source_rgba(1.0, 0.1, 0.1, 1);
//	_cairoContext->move_to(Xx,Xy);
//	_cairoContext->show_text("X");
//	_cairoContext->set_source_rgba(0.1, 1.0, 0.1, 1);
//	_cairoContext->move_to(Yx,Yy);
//	_cairoContext->show_text("Y");
//	_cairoContext->set_source_rgba(0.1, 0.1, 1.0, 1);
//	_cairoContext->move_to(Zx,Zy);
//	_cairoContext->show_text("Z");
//	_cairoContext->restore();
//	CairoSurface::syncCairoGL();
//	CairoSurface::glRender();
//      }
//
//      GLfloat Xx,Xy,Yx,Yy,Zx,Zy;
//    };


	
	/////////////////RENDER THE AXIS//////////////////////////////////////////////

	const GLdouble nearPlane = 0.1,
	  axisScale = 0.09;
    
	//The axis is in a little 100x100 pixel area in the lower left
	std::tr1::array<GLint, 4> oldviewport = context->getViewport();
	context->setViewport(0,0,100,100);
    
	GLMatrix oldproj = context->getAttachedShader()["ProjectionMatrix"].as<GLMatrix>();
	GLMatrix oldview = context->getAttachedShader()["ViewMatrix"].as<GLMatrix>();
	
	GLMatrix viewMatrix 
	  = GLMatrix::translate(0, 0, -(nearPlane + axisScale))
	  * camera.getViewRotationMatrix()
	  * GLMatrix::scale(axisScale, axisScale, axisScale);

	GLMatrix projectionMatrix
	  = GLMatrix::perspective(45, 1, nearPlane, 1000);

	context->getAttachedShader()["ViewMatrix"]  = viewMatrix;
	context->getAttachedShader()["ProjectionMatrix"] = projectionMatrix;

	context->color(0.5f,0.5f,0.5f,0.8f);
	_axis.glRender();

	context->setViewport(oldviewport);
	context->getAttachedShader()["ProjectionMatrix"] = oldproj;
	context->getAttachedShader()["ViewMatrix"] = oldview;
      }    
  }

  void Console::glRender(const magnet::GL::Camera& camera, RenderMode mode)
  {
    if (_showGrid->get_active())
      {
	using namespace magnet::GL;
	const Context::ContextPtr& context = _axis.getContext();

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

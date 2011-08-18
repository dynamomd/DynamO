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
    _cairoOverlay.init();
    initGTK();
  }

  void 
  Console::interfaceRender(const magnet::GL::Camera& camera)
  {
    //Only draw if the console has something in it or if it's visible
    if (_consoleEntries.empty() || !_visible) return;

    //Disable anything that might affect the rastering 
    glDisable(GL_DEPTH_TEST);

    using namespace magnet::GL;
    Context& context = _axis.getContext();
    //Draw the console in orthograpic projection
    context.cleanupAttributeArrays();

//    if (_showConsole->get_active())
//      {
//	float lineHeight = _consoleFont->FaceSize() / (0.5f * _viewPort->getHeight());
//	float consoleHeight = 1.0f - lineHeight;
//	
//	//Calculate how long since the last redraw
//	int tdelta = glutGet(GLUT_ELAPSED_TIME) - _glutLastTime;
//	_glutLastTime = glutGet(GLUT_ELAPSED_TIME);
//	
//	glColor3f(_consoleTextColor[0], _consoleTextColor[1], 
//		  _consoleTextColor[2]);
//	
//	glRasterPos3f(-1.0, consoleHeight, 0);
//	_consoleLayout->Render(_consoleEntries.front().second.c_str());
//	consoleHeight -= lineHeight;
//	
//	for (std::list<consoleEntry>::iterator iPtr = ++_consoleEntries.begin();
//	     iPtr != _consoleEntries.end();)
//	  {
//	    //Fade the color based on it's time in the queue
//	    glColor4f(_consoleTextColor[0], _consoleTextColor[1], 
//		      _consoleTextColor[2], 1.0f - iPtr->first / 1000.0f);
//	    glRasterPos3f(-1, consoleHeight, 0);
//	    _consoleLayout->Render(iPtr->second.c_str());
//	    iPtr->first += tdelta;
//	    consoleHeight -= lineHeight;
//	    
//	    std::list<consoleEntry>::iterator prev = iPtr++;
//	    //If this element is invisible, erase it
//	    if (prev->first > 1000) _consoleEntries.erase(prev);
//	  }
//      }

    if (_showAxis->get_active())
      {
	/////////////////RENDER THE AXIS//////////////////////////////////////////////

	const GLdouble nearPlane = 0.1,
	  axisScale = 0.09;
    
	//The axis is in a little 100x100 pixel area in the lower left
	std::tr1::array<GLint, 4> oldviewport = context.getViewport();
	context.setViewport(0,0,100,100);
    
	std::tr1::array<GLfloat, 16> oldproj 
	  = context.getAttachedShader()["ProjectionMatrix"].as<std::tr1::array<GLfloat, 16> >();
	std::tr1::array<GLfloat, 16> oldview 
	  = context.getAttachedShader()["ViewMatrix"].as<std::tr1::array<GLfloat, 16> >();
	
	GLMatrix viewMatrix 
	  = GLMatrix::translate(0, 0, -(nearPlane + axisScale))
	  * GLMatrix::rotate(camera.getTilt(), Vector(1, 0, 0))
	  * GLMatrix::rotate(camera.getPan(), Vector(0, 1, 0))
	  * GLMatrix::scale(axisScale, axisScale, axisScale);

	GLMatrix projectionMatrix
	  = GLMatrix::perspective(45, 1, nearPlane, 1000);

	context.getAttachedShader()["ViewMatrix"]  = viewMatrix;
	context.getAttachedShader()["ProjectionMatrix"] = projectionMatrix;

	context.color(0.5f,0.5f,0.5f,0.8f);
	_axis.glRender();

	_cairoOverlay.glRender(projectionMatrix * viewMatrix);

	context.setViewport(oldviewport);
	context.getAttachedShader()["ProjectionMatrix"] = oldproj;
	context.getAttachedShader()["ViewMatrix"] = oldview;
      }    

    //Restore GL state
    glEnable(GL_DEPTH_TEST);
  }

  void Console::glRender(magnet::GL::FBO&, const magnet::GL::Camera& camera, RenderMode mode)
  {

    if (_showGrid->get_active())
      {
	using namespace magnet::GL;
	Context& context = _axis.getContext();

	GLMatrix old_model_view
	  = context.getAttachedShader()["ViewMatrix"].as<std::tr1::array<GLfloat, 16> >();

	context.getAttachedShader()["ViewMatrix"]
	  = old_model_view
	  * GLMatrix::translate(camera.getViewPlanePosition())
	  * GLMatrix::rotate(-camera.getPan(), Vector(0, 1, 0))
	  * GLMatrix::rotate(-camera.getTilt(), Vector(1, 0, 0));

	context.color(1,1,1,1);
	//Back face
	context.setAttribute(Context::instanceOriginAttrIndex, 0, 0, -camera.getScreenPlaneWidth(), 0);
	context.setAttribute(Context::instanceScaleAttrIndex,
			     camera.getScreenPlaneWidth(),
			     camera.getScreenPlaneHeight(), 1);
	_grid.glRender();

	//Sides
	context.setAttribute(Context::instanceOriginAttrIndex, 
			     0.5 * camera.getScreenPlaneWidth(), 0, 
			     -0.5 * camera.getScreenPlaneWidth(), 0);
	context.rotation(M_PI / 2, Vector(0, 1, 0));
  	_grid.glRender(); //Right side
	context.setAttribute(Context::instanceOriginAttrIndex, 
			     -0.5 * camera.getScreenPlaneWidth(), 0, 
			     -0.5 * camera.getScreenPlaneWidth(), 0);
  	_grid.glRender(); //Left side

	//Top and bottom
	context.rotation(M_PI / 2, Vector(1, 0, 0));
	context.setAttribute(Context::instanceScaleAttrIndex,
			     camera.getScreenPlaneWidth(),
			     camera.getScreenPlaneWidth(), 1);
	context.setAttribute(Context::instanceOriginAttrIndex, 0,
			     -0.5 * camera.getScreenPlaneHeight(), 
			     -0.5 * camera.getScreenPlaneWidth(), 0);
	_grid.glRender();//bottom
	context.setAttribute(Context::instanceOriginAttrIndex, 0,
			     0.5 * camera.getScreenPlaneHeight(), 
			     -0.5 * camera.getScreenPlaneWidth(), 0);
	_grid.glRender();//top
	context.getAttachedShader()["ViewMatrix"] = old_model_view;
      }
  }


  void
  Console::initGTK()
  {
    _optList.reset(new Gtk::VBox);//The Vbox of options   

    {
      _showGrid.reset(new Gtk::CheckButton("Show viewing grid"));
      _showGrid->set_active(false);
      _optList->add(*_showGrid); _showGrid->show();
    }

    {
      _showConsole.reset(new Gtk::CheckButton("Show console"));
      _showConsole->set_active(false);
      _showConsole->set_sensitive(false);
      _optList->add(*_showConsole); _showConsole->show();
    }

    {
      _showAxis.reset(new Gtk::CheckButton("Show axis"));
      _showAxis->set_active(true);
      _optList->add(*_showAxis); _showAxis->show();
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

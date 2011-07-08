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

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>
#include <GL/freeglut.h>

extern const unsigned char _binary_coilfont_ttf_start[];
extern const unsigned char _binary_coilfont_ttf_end[];

namespace coil {
  void 
  Console::initOpenGL() 
  {
    _consoleFont.reset(new FTGLPixmapFont(_binary_coilfont_ttf_start,
					  _binary_coilfont_ttf_end
					  -_binary_coilfont_ttf_start));
    _consoleLayout.reset(new FTSimpleLayout());

    if (_consoleFont->Error()) 
      M_throw() << "Could not load coil's embedded font! Errno " 
		<< _consoleFont->Error();
    
    _consoleFont->FaceSize(16);

    _consoleLayout->SetFont(&(*_consoleFont));    

    if (_consoleLayout->Error()) 
      M_throw() << "Could set the font of the layout " 
		<< _consoleLayout->Error();

    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);

    resize(_viewPort->getWidth(), _viewPort->getHeight());

    _axis.init();
    _grid.init(10,10);
  }

  void 
  Console::resize(size_t width, size_t height)
  {
    _consoleLayout->SetLineLength(width);
  }

  void 
  Console::interfaceRender()
  {
    //Only draw if the console has something in it or if it's visible
    if (_consoleEntries.empty() || !_visible) return;

    //Disable anything that might affect the rastering 
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    //Draw the console in orthograpic projection
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    if (_showConsole->get_active())
      {
	float lineHeight = _consoleFont->FaceSize() / (0.5f * _viewPort->getHeight());
	float consoleHeight = 1.0f - lineHeight;
	
	//Calculate how long since the last redraw
	int tdelta = glutGet(GLUT_ELAPSED_TIME) - _glutLastTime;
	_glutLastTime = glutGet(GLUT_ELAPSED_TIME);
	
	glColor3f(_consoleTextColor[0], _consoleTextColor[1], 
		  _consoleTextColor[2]);
	
	glRasterPos3f(-1.0, consoleHeight, 0);
	_consoleLayout->Render(_consoleEntries.front().second.c_str());
	consoleHeight -= lineHeight;
	
	for (std::list<consoleEntry>::iterator iPtr = ++_consoleEntries.begin();
	     iPtr != _consoleEntries.end();)
	  {
	    //Fade the color based on it's time in the queue
	    glColor4f(_consoleTextColor[0], _consoleTextColor[1], 
		      _consoleTextColor[2], 1.0f - iPtr->first / 1000.0f);
	    glRasterPos3f(-1, consoleHeight, 0);
	    _consoleLayout->Render(iPtr->second.c_str());
	    iPtr->first += tdelta;
	    consoleHeight -= lineHeight;
	    
	    std::list<consoleEntry>::iterator prev = iPtr++;
	    //If this element is invisible, erase it
	    if (prev->first > 1000) _consoleEntries.erase(prev);
	  }
      }

    if (_showAxis->get_active())
      {
	/////////////////RENDER THE AXIS//////////////////////////////////////////////

	const GLdouble nearPlane = 0.1,
	  axisScale = 0.07;
    
	//The axis is in a little 100x100 pixel area in the lower left
	GLint viewportDim[4];
	glGetIntegerv(GL_VIEWPORT, viewportDim);
	glViewport(0,0,100,100);
    
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, 1, nearPlane, 1000.0f);
    
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    
	//near plane is at 0.1, the axis are axisScale long so
	glTranslatef (0, 0, -(nearPlane + axisScale));
    
	glColor4f (0.5f,0.5f,0.5f,0.8f); // Color the axis box a transparent white
	glBegin(GL_QUADS);		
	glVertex3f(-1,-1, 0);
	glVertex3f( 1,-1, 0);
	glVertex3f( 1, 1, 0);
	glVertex3f(-1, 1, 0);
	glEnd();
    
	glRotatef(_viewPort->getTilt(), 1.0, 0.0, 0.0);
	glRotatef(_viewPort->getPan(), 0.0, 1.0, 0.0);
	glScalef (axisScale, axisScale, axisScale);
    
	glLineWidth(2.0f);

	_axis.glRender();
    
	//Do the axis labels
	glColor3f(1,1,1);
	glRasterPos3f( 0.5,-0.5,-0.5);
	_consoleFont->Render("X");
	glRasterPos3f(-0.5, 0.5,-0.5);
	_consoleFont->Render("Y");
	glRasterPos3f(-0.5,-0.5, 0.5);
	_consoleFont->Render("Z");
	glViewport(viewportDim[0], viewportDim[1], viewportDim[2], viewportDim[3]);
      }    

    //Restore GL state
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix ();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix ();
  }

  void Console::glRender()
  {
    if (_showGrid->get_active())
      {
	//Draw a frustrum cage around the display
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	_viewPort->applyInverseModelview();

	GLfloat modelview[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

	glColor4f(1, 1, 1, 1);

	//Back face
	glLoadIdentity();
	glMultMatrixf(modelview);
	glScalef(_viewPort->getScreenPlaneWidth(), _viewPort->getScreenPlaneHeight(), 1);
	glTranslatef(0, 0, -1.0f);
	_grid.glRender();

	//Sides
	glLoadIdentity();
	glMultMatrixf(modelview);
	glScalef(_viewPort->getScreenPlaneWidth(), _viewPort->getScreenPlaneHeight(), 1);
	glRotatef(90, 0, 1, 0);
	//left
	glTranslatef(0.5f, 0, -0.5f);
	_grid.glRender();
	//right
	glTranslatef(0, 0, 1);
	_grid.glRender();

	//Top and bottom
	glLoadIdentity();
	glMultMatrixf(modelview);
	glScalef(_viewPort->getScreenPlaneWidth(), _viewPort->getScreenPlaneWidth(), 1);
	glRotatef(90, 1, 0, 0);
	//top
	glTranslatef(0, -0.5f, -0.5f *  _viewPort->getScreenPlaneHeight() / _viewPort->getScreenPlaneWidth());
	_grid.glRender();
	//bottom
	glTranslatef(0, 0, _viewPort->getScreenPlaneHeight() / _viewPort->getScreenPlaneWidth());
	_grid.glRender();

    
	glPopMatrix();
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
      _showConsole->set_active(true);
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

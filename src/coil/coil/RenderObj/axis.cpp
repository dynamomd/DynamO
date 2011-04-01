/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include <coil/RenderObj/axis.hpp>
#include <magnet/exception.hpp>

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>
#include <GL/freeglut.h>

#include <coil/glprimatives/arrow.hpp>
#include <magnet/GL/viewPort.hpp>

extern const unsigned char _binary_coilfont_ttf_start[];
extern const unsigned char _binary_coilfont_ttf_end[];

namespace coil {
  Axis::Axis():
    RenderObj("Axis")
  {}

  void 
  Axis::initOpenGL() 
  {
    _axisFont.reset(new FTGLPixmapFont(_binary_coilfont_ttf_start,
				       _binary_coilfont_ttf_end
				       -_binary_coilfont_ttf_start));
    if (_axisFont->Error()) 
      M_throw() << "Could not load coil's embedded font! Errno " 
		<< _axisFont->Error();

  }

  void 
  Axis::interfaceRender()
  {
    GLdouble nearPlane = 0.1,
      axisScale = 0.07;
   
    //Only draw if the console has something in it or if it's visible
    if (!_visible) return;

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    //We want the arrow drawing to always succeed
    glDisable(GL_DEPTH_TEST);
    
    //The axis is in a little 100x100 pixel area in the lower left
    GLint viewportDim[4];
    glGetIntegerv(GL_VIEWPORT, viewportDim);
    glViewport(0,0,100,100);
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(45.0f, 1, nearPlane, 1000.0f);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity();
    
    //near plane is at 0.1, the axis are axisScale long so
    glTranslatef (0, 0, -(nearPlane + axisScale));
    
    glColor4f (4.0/256,104.0/256.0,202.0/256.0, 0.5); // Color the axis box a transparent blue
    glBegin(GL_QUADS);		
    glVertex3f(-1,-1, 0);
    glVertex3f( 1,-1, 0);
    glVertex3f( 1, 1, 0);
    glVertex3f(-1, 1, 0);
    glEnd();
    
    glRotatef(_viewPort->getTilt(), 1.0, 0.0, 0.0);
    glRotatef(_viewPort->getPan(), 0.0, 1.0, 0.0);
    glScalef (axisScale, axisScale, axisScale);
    
    glLineWidth (2.0);
    
    glColor3f (1,0,0); // X axis is red.
    coil::glprimatives::drawArrow(Vector(-0.5,-0.5,-0.5),
				  Vector( 0.5,-0.5,-0.5));
    
    glColor3f (0,1,0); // Y axis is green.
    coil::glprimatives::drawArrow(Vector(-0.5,-0.5,-0.5), 
				  Vector(-0.5, 0.5,-0.5));
    
    glColor3f (0,0,1); // Z axis is blue.
    coil::glprimatives::drawArrow(Vector(-0.5,-0.5,-0.5),
				  Vector(-0.5,-0.5, 0.5));
    
    //Do the axis labels
    glColor3f(1,1,1);
    _axisFont->FaceSize(16);
    
    glRasterPos3f( 0.5,-0.5,-0.5);
    _axisFont->Render("X");
    glRasterPos3f(-0.5, 0.5,-0.5);
    _axisFont->Render("Y");
    glRasterPos3f(-0.5,-0.5, 0.5);
    _axisFont->Render("Z");
    
    //Restore GL state
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();    
    glViewport(viewportDim[0], viewportDim[1], viewportDim[2], viewportDim[3]);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
  }
  
}

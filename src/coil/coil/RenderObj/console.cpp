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

#include <coil/RenderObj/console.hpp>
#include <magnet/exception.hpp>
#include <magnet/clamp.hpp>

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>
#include <GL/freeglut.h>

extern const unsigned char _binary_src_coil_coil_coilfont_ttf_start[];
extern const unsigned char _binary_src_coil_coil_coilfont_ttf_end[];

namespace coil {
  Console::Console(size_t width, size_t height, float r, float g, float b):
    RenderObj("Console"),
    _width(width),
    _height(height)
  {
    _color[0] = r;
    _color[1] = g;
    _color[2] = b;
  }

  void 
  Console::initOpenGL() 
  {
    _consoleFont.reset(new FTGLPixmapFont(_binary_src_coil_coil_coilfont_ttf_start,
					  _binary_src_coil_coil_coilfont_ttf_end
					  -_binary_src_coil_coil_coilfont_ttf_start));
    _consoleLayout.reset(new FTSimpleLayout());

    if (_consoleFont->Error()) 
      M_throw() << "Could not load coil's embedded font! Errno " 
		<< _consoleFont->Error();
    
    _consoleLayout->SetFont(&(*_consoleFont));
    
    if (_consoleLayout->Error()) 
      M_throw() << "Could set the font of the layout " 
		<< _consoleLayout->Error();

    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);

    resize(_width, _height);
  }

  void 
  Console::resize(size_t width, size_t height)
  {
    _width = width;
    _height = height;
    _consoleLayout->SetLineLength(_width);
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
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    float lineHeight = _consoleFont->FaceSize() / (0.5f * _height);
    float consoleHeight = 1.0f - lineHeight;

    //Calculate how long since the last redraw
    int tdelta = glutGet(GLUT_ELAPSED_TIME) - _glutLastTime;
    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);

    glColor4f(_color[0], _color[1], _color[2], 1);
    glRasterPos3f(-1.0, consoleHeight, 0);
    _consoleLayout->Render(_consoleEntries.front().second.c_str());
    consoleHeight -= lineHeight;
  
    for (std::list<consoleEntry>::iterator iPtr = ++_consoleEntries.begin();
	 iPtr != _consoleEntries.end();)
      {
	//Fade the color based on it's time in the queue
	glColor4f(_color[0], _color[1], _color[2],  1.0f - iPtr->first / 1000.0f);
	glRasterPos3f(-1, consoleHeight, 0);
	_consoleLayout->Render(iPtr->second.c_str());
	iPtr->first += tdelta;
	consoleHeight -= lineHeight;

	std::list<consoleEntry>::iterator prev = iPtr++;
	//If this element is invisible, erase it
	if (prev->first > 1000) _consoleEntries.erase(prev);
      }

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
  }
  
}

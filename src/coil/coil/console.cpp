/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <coil/console.hpp>
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
  Console::Console(size_t width, size_t height):
    _consoleFont(_binary_src_coil_coil_coilfont_ttf_start,
		 _binary_src_coil_coil_coilfont_ttf_end
		 -_binary_src_coil_coil_coilfont_ttf_start)
  {
    if (_consoleFont.Error()) 
      M_throw() << "Could not load coil's embedded font! Errno " 
		<< _consoleFont.Error();

    _consoleFont.FaceSize(12);

    _consoleLayout.SetFont(&_consoleFont);

    if (_consoleLayout.Error()) 
      M_throw() << "Could set the font of the layout " 
		<< _consoleLayout.Error();

    _glutLastTime = glutGet(GLUT_ELAPSED_TIME);

    resize(width, height);
  }

  void 
  Console::resize(size_t width, size_t height)
  {
    _width = width;
    _height = height;
    _consoleLayout.SetLineLength(_width);
  }

  void 
  Console::draw()
  {
    //Disable anything that might affect the rastering 
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    //Draw the console in orthograpic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    GLfloat lineHeight = _consoleFont.FaceSize() / (0.5f * _height);
    GLfloat consoleHeight = 1.0f - lineHeight;
    //Calculate how long since the last redraw
    int glutTime = glutGet(GLUT_ELAPSED_TIME);
    int tdelta = glutTime - _glutLastTime;
    _glutLastTime = glutTime;
  
    //glEnable(GL_BLEND);
    for (std::list<consoleEntry>::iterator iPtr = _consoleEntries.begin();
	 iPtr != _consoleEntries.end(); ++iPtr)
      {
	std::ostringstream os;
	GLfloat textAlpha = iPtr->first / 1000;
	os << iPtr->second << " Alpha val " << textAlpha;
	std::string str = os.str();
	glColor4f(0, 0, 0, textAlpha);
	glRasterPos3f(-1.0, consoleHeight, 0);
	_consoleLayout.Render(str.c_str());
	iPtr->first += tdelta;
	consoleHeight -= lineHeight;
      }

    //glDisable(GL_BLEND);
  }
  
}

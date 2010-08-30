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

#pragma once
#include <sstream>
#include <GL/freeglut.h>

#define GL_FONT_FACTORY(F)\
  F(BITMAP_8_BY_13, GLUT_BITMAP_8_BY_13)			\
  F(BITMAP_9_BY_15, GLUT_BITMAP_9_BY_15)			\
  F(BITMAP_TIMES_ROMAN_10, GLUT_BITMAP_TIMES_ROMAN_10)		\
  F(BITMAP_TIMES_ROMAN_24, GLUT_BITMAP_TIMES_ROMAN_24)		\
  F(BITMAP_HELVETICA_10, GLUT_BITMAP_HELVETICA_10)		\
  F(BITMAP_HELVETICA_12, GLUT_BITMAP_HELVETICA_12)		\
  F(BITMAP_HELVETICA_18, GLUT_BITMAP_HELVETICA_18)

class GLScribe 
{
public:
  static GLScribe cout;
  
  //Font enum building
#define buildEnum(NAME,VAL) \
  NAME,

  typedef enum {
    GL_FONT_FACTORY(buildEnum)
  } font; 
#undef buildEnum

  //A manipulator to change the plot position
  struct cursor
  {
    cursor(GLfloat x, GLfloat y, GLfloat z):
      _x(x),_y(y),_z(z)
    {}
    
    GLfloat _x, _y, _z;
  };

  template<class T>
  GLScribe& operator<<(const T& data)
  {
    //We reuse the _stream member as it stores the state of the stream (setprecision etc.)
    //Get the string representation
    _stream << data;
    std::string output = _stream.str();
    //Blank the output stream
    _stream.str("");
    
    //Print it on the screen
    return operator<<<std::string>(output);
  }

  void* _font;
    
  private:
  GLScribe():
    _font(GLUT_BITMAP_8_BY_13)
  {}
  
  std::ostringstream _stream;
};

template<>
inline GLScribe& 
GLScribe::operator<<<std::string>(const std::string& data)
{
  for (size_t i = 0; i < data.size(); ++i)
    glutBitmapCharacter(_font, data[i]);

  return *this;
}

template<>
inline GLScribe& 
GLScribe::operator<<<GLScribe::font>(const GLScribe::font& val)
{
#define buildSwitch(NAME,VAL) \
  case NAME: _font = VAL; break;

  switch (val)
    {
      GL_FONT_FACTORY(buildSwitch)
    }

#undef buildSwitch

  return *this;
}

template<>
inline GLScribe& 
GLScribe::operator<<<GLScribe::cursor>(const GLScribe::cursor& newpos)
{
  glRasterPos3f(newpos._x, newpos._y, newpos._z);

  return *this;
}

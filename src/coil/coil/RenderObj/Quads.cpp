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
#include "Quads.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

RQuads::RQuads(std::string name):
  RTriangles(name) {}

void 
RQuads::glRender()
{
  if (!_visible) return;

  if (!_colBuff.empty())
    {
      _colBuff.bind(magnet::GL::ARRAY);
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
      glEnableClientState(GL_COLOR_ARRAY);
    }

  if (!_normBuff.empty())
    {
      _normBuff.bind(magnet::GL::ARRAY);
      glNormalPointer(GL_FLOAT, 0, 0);
      glEnableClientState(GL_NORMAL_ARRAY); 
    }
  
  _posBuff.bind(magnet::GL::ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  _elementBuff.bind(magnet::GL::ELEMENT_ARRAY);
  
  glEnableClientState(GL_VERTEX_ARRAY);
  
  switch (_RenderMode)
    {
    case TRIANGLES:
      glDrawElements(GL_QUADS, _elementBuff.size(), GL_UNSIGNED_INT, 0);
      break;
    case LINES:
      glDrawElements(GL_LINES, _elementBuff.size(), GL_UNSIGNED_INT, 0);
      break;
    case POINTS:
      glDrawElements(GL_POINTS, _elementBuff.size(), GL_UNSIGNED_INT, 0);
      break;
    }
 
  if (_colBuff.size())
    glDisableClientState(GL_COLOR_ARRAY);	
  if (_normBuff.size())
    glDisableClientState(GL_NORMAL_ARRAY);

  glDisableClientState(GL_VERTEX_ARRAY);

  if (_renderNormals && _normBuff.size())
    {
      const GLfloat* posPointer = _posBuff.map();
      const GLfloat* normPointer = _normBuff.map();

      const float scale = 0.005;
      for (size_t i= 0; i < _posBuff.size(); i+= 3)
	{
	  Vector point1, point2;
	  for (size_t iDim = 0; iDim < 3; ++iDim)
	    {
	      point1[iDim] = posPointer[i + iDim];
	      point2[iDim] = point1[iDim] + scale * normPointer[i + iDim];
	    }
	  coil::glprimatives::drawArrow(point1, point2);
	}
      
      _posBuff.unmap();
      _normBuff.unmap();
    }

}

void 
RQuads::setGLElements(std::vector<GLuint>& Elements)
{
  if (!Elements.size())
    throw std::runtime_error("Elements.size() == 0!");

  if (Elements.size() % 4) 
    throw std::runtime_error("Elements.size() is not a multiple of 4!");

  _elementBuff.init(Elements);
}

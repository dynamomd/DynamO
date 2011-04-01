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
#include "Quads.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

RQuads::RQuads(std::string name):
  RTriangles(name) {}

void 
RQuads::glRender()
{
  if (!_visible) return;
  glBindBufferARB(GL_ARRAY_BUFFER, _colBuff);
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);

  glBindBufferARB(GL_ARRAY_BUFFER, _posBuff);
  glVertexPointer(3, GL_FLOAT, 0, 0);
  
  glBindBufferARB(GL_ARRAY_BUFFER, _normBuff);
  glNormalPointer(GL_FLOAT, 0, 0);
  
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, _elementBuff);
  
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY); 
  glEnableClientState(GL_VERTEX_ARRAY);
  
  switch (_RenderMode)
    {
    case TRIANGLES:
      glDrawElements(GL_QUADS, _elementBuffSize, GL_UNSIGNED_INT, 0);
      break;
    case LINES:
      glDrawElements(GL_LINES, _elementBuffSize, GL_UNSIGNED_INT, 0);
      break;
    case POINTS:
      glDrawElements(GL_POINTS, _elementBuffSize, GL_UNSIGNED_INT, 0);
      break;
    }
 
  glDisableClientState(GL_COLOR_ARRAY);	
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  if (_renderNormals)
    {
      glBindBuffer(GL_ARRAY_BUFFER, _posBuff);
      const float* posPointer = (const float*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
      glBindBuffer(GL_ARRAY_BUFFER, _normBuff);
      const float* normPointer = (const float*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

      const float scale = 0.005;
      for (size_t i= 0; i < _posBuffSize; i+= 3)
	{
	  Vector point1, point2;
	  for (size_t iDim = 0; iDim < 3; ++iDim)
	    {
	      point1[iDim] = posPointer[i + iDim];
	      point2[iDim] = point1[iDim] + scale * normPointer[i + iDim];
	    }
	  coil::glprimatives::drawArrow(point1, point2);
	}
      
      glBindBuffer(GL_ARRAY_BUFFER, _posBuff);
      glUnmapBuffer(GL_ARRAY_BUFFER);
      glBindBuffer(GL_ARRAY_BUFFER, _normBuff);
      glUnmapBuffer(GL_ARRAY_BUFFER);
    }

}

void 
RQuads::setGLElements(std::vector<int>& Elements)
{
  if (!Elements.size())
    throw std::runtime_error("Elements.size() == 0!");

  if (Elements.size() % 4) 
    throw std::runtime_error("Elements.size() is not a multiple of 4!");

  if (_elementBuffSize)
    glDeleteBuffers(1, &_elementBuff);
  
  _elementBuffSize = Elements.size();

  glGenBuffersARB(1, &_elementBuff);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, _elementBuff);
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, Elements.size() * sizeof(int), &Elements[0], 
	       GL_STATIC_DRAW);
}

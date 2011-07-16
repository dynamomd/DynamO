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
    _colBuff.attachToColor(4);

  if (_normBuff.size())
    _normBuff.attachToNormal();
  
  _posBuff.getContext().cleanupAttributeArrays();
  _posBuff.attachToVertex(3);
  
  switch (_RenderMode)
    {
    case TRIANGLES:
      _elementBuff.drawElements(magnet::GL::element_type::QUADS);
      break;
    case LINES:
      _elementBuff.drawElements(magnet::GL::element_type::LINES);
      break;
    case POINTS:
      _elementBuff.drawElements(magnet::GL::element_type::POINTS);
      break;
    }
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

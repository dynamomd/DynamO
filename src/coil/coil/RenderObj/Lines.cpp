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
#include "Lines.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

RLines::RLines(size_t N, std::string name):
  RenderObj(name), _N(N)
{}

RLines::~RLines()
{
  releaseCLGLResources();
}

void 
RLines::initOpenGL()
{
  {//Setup initial vertex positions
    std::vector<float> VertexPos(3 * _N * 2, 0.0);
    for (size_t i(0); i < _N; ++i)
      { 
	VertexPos[6*i+0] = i * 1.0f / _N;
	VertexPos[6*i+1] = i * 1.0f / _N;
	VertexPos[6*i+2] = i * 1.0f / _N;
	VertexPos[6*i+3] = i * 1.0f / _N;
	VertexPos[6*i+4] = (i + 0.5f) * 1.0f / _N;
	VertexPos[6*i+5] = i * 1.0f / _N;
      }
    setGLPositions(VertexPos);
  }
  
  {
    std::vector<cl_uchar4> VertexColor(_N * 2);
    
    for (size_t icol = 0; icol < _N * 2; ++icol)
      {
	VertexColor[icol].s[0] = 255;
	VertexColor[icol].s[1] = 255;
	VertexColor[icol].s[2] = 255;
	VertexColor[icol].s[3] = 255;
      }

    setGLColors(VertexColor);
  }
   
  {//Setup initial element data
    std::vector<GLuint> ElementData(2 * _N, 0);

    for (size_t i(0); i < _N; ++i)
      {
	ElementData[2 * i + 0] = 2 * i + 0;
	ElementData[2 * i + 1] = 2 * i + 1;
      }
    
    setGLElements(ElementData);
  }
}

void 
RLines::initOpenCL()
{
  initOCLVertexBuffer(magnet::GL::Context::getContext().getCLContext());
}

void 
RLines::glRender()
{
  if (!_visible) return;
  
  _colBuff.bind(magnet::GL::ARRAY);
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
  
  _posBuff.bind(magnet::GL::ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);
  
  _elementBuff.bind(magnet::GL::ELEMENT_ARRAY);
  
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  
  switch (_RenderMode)
    {
    case TRIANGLES:
    case LINES:
      glDrawElements(GL_LINES, _elementBuff.size(), GL_UNSIGNED_INT, 0);
      break;
    case POINTS:
      glDrawElements(GL_POINTS, _elementBuff.size(), GL_UNSIGNED_INT, 0);
      break;
    }
 
  glDisableClientState(GL_COLOR_ARRAY);	
  glDisableClientState(GL_VERTEX_ARRAY);
}

void 
RLines::setGLColors(std::vector<cl_uchar4>& VertexColor)
{
  if (!VertexColor.size())
    throw std::runtime_error("VertexColor.size() == 0!");

  if (!_posBuff.empty())
    if ((VertexColor.size()) != (_posBuff.size() / 3))
      throw std::runtime_error("VertexColor.size() != posBuffSize/3");

  _colBuff.init(VertexColor, magnet::GL::STREAM_DRAW);
}

void 
RLines::setGLPositions(std::vector<float>& VertexPos)
{
  if (!VertexPos.size())
    throw std::runtime_error("VertexPos.size() == 0!");
  
  if (VertexPos.size() % 3)
    throw std::runtime_error("VertexPos.size() is not a multiple of 3!");

  if (!_colBuff.empty())
    if ((_colBuff.size()) != (VertexPos.size() / 3))
      throw std::runtime_error("VertexPos.size()/3 != colBuffSize/4 ");
  
  _posBuff.init(VertexPos, magnet::GL::STREAM_DRAW);
}

void 
RLines::initOCLVertexBuffer(cl::Context Context)
{
  _clbuf_Positions = cl::GLBuffer<GLfloat>(Context, CL_MEM_READ_WRITE, _posBuff);
}

void 
RLines::initOCLColorBuffer(cl::Context Context)
{
  _clbuf_Colors = cl::GLBuffer<cl_uchar4>(Context, CL_MEM_READ_WRITE, _colBuff);
}

void 
RLines::initOCLElementBuffer(cl::Context Context)
{
  _clbuf_Elements = cl::GLBuffer<GLuint>(Context, CL_MEM_READ_WRITE, _elementBuff);
}

void 
RLines::setGLElements(std::vector<GLuint>& Elements)
{
  if (!Elements.size())
    throw std::runtime_error("Elements.size() == 0!");

  if (Elements.size() % 2)
    throw std::runtime_error("Elements.size() is not a multiple of 3!");

  _elementBuff.init(Elements);
}
 
void 
RLines::releaseCLGLResources()
{
  _colBuff.deinit();
  _posBuff.deinit();
  _elementBuff.deinit();
}

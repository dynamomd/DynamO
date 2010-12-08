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
#include "Lines.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>

RLines::RLines(size_t N):
  _N(N),
  _colBuffSize(0),
  _posBuffSize(0),
  _elementBuffSize(0)
{}

RLines::~RLines()
{
  if (_colBuffSize)
    glDeleteBuffersARB(1, &_colBuff);

  if (_posBuffSize)
    glDeleteBuffersARB(1, &_posBuff);

  if (_elementBuffSize)
    glDeleteBuffersARB(1, &_elementBuff);
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
    std::vector<int> ElementData(2 * _N, 0);

    for (size_t i(0); i < _N; ++i)
      {
	ElementData[2 * i + 0] = 2 * i + 0;
	ElementData[2 * i + 1] = 2 * i + 1;
      }
    
    setGLElements(ElementData);
  }

}

void 
RLines::initOpenCL(magnet::CL::CLGLState& CLState)
{
  initOCLVertexBuffer(CLState.getContext());
}

void 
RLines::glRender()
{
  if (!_visible) return;

  glBindBufferARB(GL_ARRAY_BUFFER, _colBuff);
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);

  glBindBufferARB(GL_ARRAY_BUFFER, _posBuff);
  glVertexPointer(3, GL_FLOAT, 0, 0);
  
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, _elementBuff);
  
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  
  switch (_RenderMode)
    {
    case TRIANGLES:
    case LINES:
      glDrawElements(GL_LINES, _elementBuffSize, GL_UNSIGNED_INT, 0);
      break;
    case POINTS:
      glDrawElements(GL_POINTS, _elementBuffSize, GL_UNSIGNED_INT, 0);
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

  if (_posBuffSize)
    if ((VertexColor.size()) != (_posBuffSize / 3))
      throw std::runtime_error("VertexColor.size() != posBuffSize/3");

  if (_colBuffSize)
    glDeleteBuffers(1, &_colBuff);

  _colBuffSize = VertexColor.size();
  
  glGenBuffersARB(1, &_colBuff);
  
  glBindBufferARB(GL_ARRAY_BUFFER, _colBuff);
  glBufferDataARB(GL_ARRAY_BUFFER, VertexColor.size() * sizeof(cl_uchar4), &VertexColor[0], 
	       GL_STREAM_DRAW);
}

void 
RLines::setGLPositions(std::vector<float>& VertexPos)
{
  if (!VertexPos.size())
    throw std::runtime_error("VertexPos.size() == 0!");
  
  if (VertexPos.size() % 3)
    throw std::runtime_error("VertexPos.size() is not a multiple of 3!");

  if (_colBuffSize)
    if ((_colBuffSize) != (VertexPos.size() / 3))
      throw std::runtime_error("VertexPos.size()/3 != colBuffSize/4 ");
  
  if (_posBuffSize)
    glDeleteBuffers(1, &_posBuff);

  _posBuffSize = VertexPos.size();

  glGenBuffersARB(1, &_posBuff);

  glBindBufferARB(GL_ARRAY_BUFFER, _posBuff);
  glBufferDataARB(GL_ARRAY_BUFFER, VertexPos.size() * sizeof(float), &VertexPos[0], 
	       GL_STREAM_DRAW);
}

void 
RLines::initOCLVertexBuffer(cl::Context Context)
{
  _clbuf_Positions = cl::GLBuffer(Context, CL_MEM_READ_WRITE, _posBuff, GL_ARRAY_BUFFER);
}

void 
RLines::initOCLColorBuffer(cl::Context Context)
{
  _clbuf_Colors = cl::GLBuffer(Context, CL_MEM_READ_WRITE, _colBuff, GL_ARRAY_BUFFER);
}

void 
RLines::initOCLElementBuffer(cl::Context Context)
{
  _clbuf_Elements = cl::GLBuffer(Context, CL_MEM_READ_WRITE, _elementBuff, GL_ELEMENT_ARRAY_BUFFER);
}

void 
RLines::setGLElements(std::vector<int>& Elements)
{
  if (!Elements.size())
    throw std::runtime_error("Elements.size() == 0!");

  if (Elements.size() % 2)
    throw std::runtime_error("Elements.size() is not a multiple of 3!");

  if (_elementBuffSize)
    glDeleteBuffers(1, &_elementBuff);
  
  _elementBuffSize = Elements.size();

  glGenBuffersARB(1, &_elementBuff);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, _elementBuff);
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, Elements.size() * sizeof(int), &Elements[0], 
	       GL_STATIC_DRAW);
}

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
#pragma once
#include <magnet/CL/CLGL.hpp>
#include <coil/Maths/Maths.h>
#include <coil/extcode/vector2.hpp>
#include <magnet/GL/viewPort.hpp>

class RenderObj
{
public:
  RenderObj():
    _RenderMode(TRIANGLES),
    _renderNormals(false),
    _visible(true)
  {}
  
  ~RenderObj() {}

  virtual void initOpenGL() {}
  virtual void initOpenCL(magnet::CL::CLGLState&) {}

  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&) {}
  virtual void glRender() {}
  virtual void interfaceRender() {}

  virtual void initPicking(magnet::CL::CLGLState& CLState, cl_uint& offset) {}
  virtual void pickingRender() {}
  virtual void finishPicking(magnet::CL::CLGLState& CLState, cl_uint& offset, const cl_uint val) {}

  virtual void resize(size_t width, size_t height) {}

  enum RenderModeType 
    {
      POINTS,
      LINES,
      TRIANGLES
    };

  void setRenderMode(RenderModeType rm) { _RenderMode = rm; } 
  
  inline void setDisplayNormals(bool val) { _renderNormals = val; }
  inline void setVisible(bool val) { _visible = val; }

protected:
  RenderModeType _RenderMode;

  bool _renderNormals;
  bool _visible;
};

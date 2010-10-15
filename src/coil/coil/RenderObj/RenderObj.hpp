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

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>


#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <magnet/CL/CLGL.hpp>
#include <coil/Maths/Maths.h>
#include <coil/extcode/vector2.hpp>


#include <magnet/GL/viewPort.hpp>

class RenderObj
{
public:
  RenderObj();
  ~RenderObj();

  virtual void initOpenGL() = 0;
  virtual void initOpenCL(magnet::CL::CLGLState&) = 0;

  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&) = 0;
  virtual void glRender() = 0;

  enum RenderModeType 
    {
      POINTS,
      LINES,
      TRIANGLES
    };

  void setRenderMode(RenderModeType rm) { _RenderMode = rm; } 

protected:
  RenderModeType _RenderMode;
};

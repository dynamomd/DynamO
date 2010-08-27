#pragma once

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>


#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

class RenderObj
{
public:
  RenderObj(bool hostTransfers);
  ~RenderObj();

  virtual void clTick(cl::CommandQueue& CmdQ, cl::Context& Context) = 0;
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
  bool _hostTransfers;

};

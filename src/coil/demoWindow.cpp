#define GL_GLEXT_PROTOTYPES
#include "demoWindow.hpp"
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include "GLBuffer.hpp"


demoWindow::demoWindow(GlutMaster& gMaster,
                       int setWidth, int setHeight,
                       int setInitPositionX, int setInitPositionY,
                       const char * title,
		       cl::Platform& plat
		       ):
  CLGLWindow(gMaster, setWidth, setHeight, setInitPositionX, setInitPositionY, title, plat)
{
}

void demoWindow::CallBackDisplayFunc(void)
{ 
  CLGLWindow::CallBackDisplayFunc();  
}







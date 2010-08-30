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

#define GL_GLEXT_PROTOTYPES
#include "clWindow.hpp"
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>

inline float clamp(float x, float a, float b)
{
    return x < a ? a : (x > b ? b : x);
}

//a is the head, b the tail
void drawArrow(Vector a, Vector b)
{
  Vector arrowAxis = a - b;
  Vector headpoint = b + arrowAxis * 0.75;
  Vector headaxis = (arrowAxis ^ Vector(1,0,0));
  double headaxisnorm = headaxis.nrm();
  if (headaxisnorm == 0)
    {
      headaxis = ((a-b) ^ Vector(0,0,1));
      headaxisnorm = headaxis.nrm();
    }

  headaxis *= 0.15 * arrowAxis.nrm() / headaxisnorm;

  GLfloat A[]={a.x, a.y, a.z},
    B[]={b.x, b.y, b.z},
      C[]={headpoint.x + headaxis.x, headpoint.y + headaxis.y, headpoint.z + headaxis.z},
    D[]={headpoint.x - headaxis.x, headpoint.y - headaxis.y, headpoint.z - headaxis.z};
    
  glBegin (GL_LINES);
  glVertex3fv (A);
  glVertex3fv (B);  
  glVertex3fv (A);
  glVertex3fv (C);
  glVertex3fv (A);
  glVertex3fv (D);
  glEnd();
}



#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "GLBuffer.hpp"
#include <include/glscribe.hpp>

CLGLWindow::CLGLWindow(GlutMaster& gMaster,
                       int setWidth, int setHeight,
                       int initPosX, int initPosY,
                       std::string title,
		       cl::Platform& plat,
		       bool hostTransfers
		       ):
  _clplatform(plat),
  _height(setHeight),
  _width(setWidth),
  _glutMaster(gMaster),
  keyState(DEFAULT),
  windowTitle(title),
  FPSmode(false),
  frameCounter(0),
  _rotatex(0),
  _rotatey(0),
  _cameraX(0),
  _cameraY(0),
  _cameraZ(0),
  _mouseSensitivity(0.3),
  _moveSensitivity(0.005),
  _specialKeys(0),
  _hostTransfers(hostTransfers)
{
  for (size_t i(0); i < 256; ++i) keyStates[i] = false;
  
  initOpenGL(initPosX, initPosY);
  initOpenCL();
}

CLGLWindow::~CLGLWindow()
{
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    delete *iPtr;
}

void 
CLGLWindow::CameraSetup()
{
  float moveAmp = (_currFrameTime - _lastFrameTime) * _moveSensitivity;

  int _forward = (keyStates['w'] || keyStates['W']) - (keyStates['s'] || keyStates['S']);
  int _sideways = (keyStates['d'] || keyStates['D']) - (keyStates['a'] || keyStates['A']);
  int _vertical = (keyStates['q'] || keyStates['Q']) - (keyStates['z'] || keyStates['Z']);

  //Forward/Backward movement
  _cameraZ -= _forward * moveAmp * std::cos(_rotatey * (M_PI/ 180)) 
    * std::sin(_rotatex  * (M_PI/ 180) + M_PI * 0.5);  
  _cameraX -= _forward * moveAmp * std::cos(_rotatey * (M_PI/ 180)) 
    * std::cos(_rotatex  * (M_PI/ 180) + M_PI * 0.5);
  _cameraY += -_forward * moveAmp * std::sin(_rotatey * (M_PI/ 180));

  //Strafe movement
  _cameraZ += _sideways * moveAmp * std::sin(_rotatex * (M_PI/ 180));
  _cameraX += _sideways * moveAmp * std::cos(_rotatex * (M_PI/ 180));

  //Vertical movement
  _cameraY += _vertical * moveAmp;

  glLoadIdentity();
  //gluLookAt(-viewscale, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  glRotatef(_rotatey, 1.0, 0.0, 0.0);
  glRotatef(_rotatex, 0.0, 1.0, 0.0);
  glTranslatef(-_cameraX,-_cameraY,-_cameraZ);

  GLfloat light0_diffuse[] = {1.0, 1.0, 1.0, 1.0}; //diffuse intensity of the light
  GLfloat light0_ambient[] = {0.3, 0.3, 0.3, 1.0}; 
  GLfloat light0_position[] = {0.0, 0.0, -2.0, 0.0};
  GLfloat light1_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat light1_ambient[] = {0.3, 0.3, 0.3, 1.0}; 
  GLfloat light1_position[] = {-2.0, -3.0, -1.0, 0.0};
   
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
   
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  //glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  //glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
  //glLightfv(GL_LIGHT1, GL_POSITION, light1_position);


  _cameraDirection = Rodrigues(Vector(0,-_rotatex * M_PI/180,0)) * Rodrigues(Vector(-_rotatey * M_PI/180.0,0,0)) * Vector(0,0,-1);
  
}


void 
CLGLWindow::initOpenGL(int initPosX, int initPosY)
{
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowSize(_width, _height);
  glutInitWindowPosition(initPosX, initPosY);

  _glutMaster.CallGlutCreateWindow(windowTitle.c_str(), this);
  glViewport(0, 0, _width, _height);   // This may have to be moved to after the next line on some platforms

  if (glewInit() != GLEW_OK)
    std::runtime_error("Failed initialising GLEW (GL Extension Wrangler)");

  if (!glewIsSupported("GL_VERSION_2_0 GL_ARB_pixel_buffer_object"))
    std::cout << "WARNING: ARB Pixel Buffer Objects are not supported!\n"
      "WARNING: Maybe due to indirect rendering but probably because you have a poor Graphics Card/Driver.\n"
      "WARNING: Continuing anyway as we don't manipulate pixel data, yet!";

  if (!glewIsSupported("GL_VERSION_2_0 GL_ARB_vertex_buffer_object"))
    std::runtime_error("Vertex Buffer Objects are not supported by your GPU/Driver, sorry."); 

  //Now begins the example
  glClearColor(0.8,0.8,0.8,1.0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  //Both the front and back materials track the current color
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL); //and enable it

  glEnable(GL_BLEND); //Enable blending
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Blend colors using the alpha channel

  //Light our scene!
  glEnable(GL_LIGHTING);

  //Setup the viewport
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  CallBackReshapeFunc(_width, _height);

  glMatrixMode(GL_MODELVIEW);
  CameraSetup();


  glShadeModel(GL_SMOOTH);

  glEnable(GL_CULL_FACE);//Speed rendering by culling faces (only slight speedup)

  //Setup the keyboard controls
  glutIgnoreKeyRepeat(1);

  //Finally, make this window the idle one
  _glutMaster.SetIdleToCurrentWindow();
  _glutMaster.EnableIdleFunction();

  _currFrameTime = glutGet(GLUT_ELAPSED_TIME);


}

void 
CLGLWindow::initOpenCL()
{
  //Ok Now init OpenCL.
  //Create a OpenCL context from an OpenGL one
  
  GLXContext GLContext = glXGetCurrentContext();
  
  std::cout << "Attempting to make a shared OpenCL/OpenGL context.....";
  if (GLContext == NULL)
    throw std::runtime_error("Failed to obtain the GL context");
  
  cl_context_properties cpsGL[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)_clplatform(),
				    CL_GLX_DISPLAY_KHR, (intptr_t) glXGetCurrentDisplay(),
				    CL_GL_CONTEXT_KHR, (intptr_t) GLContext, 0};
  
  //create the OpenCL context 
  try {
    _clcontext = cl::Context(CL_DEVICE_TYPE_ALL, cpsGL, NULL, NULL);
    std::cout << "Success!\n";
  } catch(cl::Error& err)
    {
      std::cout << "\nFailed to create an OpenCL context from the OpenGL one.\n"
		<< "Try a selecting a different OpenCL platform or a newer driver!\n"
		<< std::endl;
    }
  
  if (_clcontext() == NULL)
    {
      _hostTransfers = true;
      
      std::cout << "Attempting to create a standard OpenCL context. This will force host transfers on.\n";
      cl_context_properties cpsFallBack[] = {CL_CONTEXT_PLATFORM, 
					     (cl_context_properties)_clplatform(),
					     0};
      try {
	_clcontext = cl::Context(CL_DEVICE_TYPE_ALL, 
				 cpsFallBack, NULL, NULL);
      } catch (cl::Error& err)
	{
	  throw std::runtime_error("Failed to create a normal OpenCL context from the supplied platform.");
	}
    }

  if (_hostTransfers) 
    std::cout << "Host transfers have been enabled, slow performance is expected\n";

  //Grab the first device
  std::vector<cl::Device> devices = _clcontext.getInfo<CL_CONTEXT_DEVICES>();

  _cldevice = devices.front();

  std::cout << "Found these usable OpenCL Devices\n";
  for (std::vector<cl::Device>::const_iterator iPtr = devices.begin();
       iPtr != devices.end(); ++iPtr)
    switch(iPtr->getInfo<CL_DEVICE_TYPE>()) 
      {
      case CL_DEVICE_TYPE_ACCELERATOR:
	std::cout << " ACCELERATOR:" << iPtr->getInfo<CL_DEVICE_NAME>()  << "\n";
	break;
      case CL_DEVICE_TYPE_CPU:
	std::cout << " CPU:" << iPtr->getInfo<CL_DEVICE_NAME>() << "\n";
	break;
      case CL_DEVICE_TYPE_GPU:
	std::cout << " GPU:" << iPtr->getInfo<CL_DEVICE_NAME>() << "\n";
	_cldevice = *iPtr;
	break;
      default:
	std::cout << " DEFAULT:" << iPtr->getInfo<CL_DEVICE_NAME>() << "\n";
      }

  //Just check if there is a GPU to use instead of a CPU
  std::cout << "\nUsing OpenCL Device ";
  switch(_cldevice.getInfo<CL_DEVICE_TYPE>()) 
    {
    case CL_DEVICE_TYPE_ACCELERATOR:
      std::cout << " ACCELERATOR:" << _cldevice.getInfo<CL_DEVICE_NAME>();
      break;
    case CL_DEVICE_TYPE_CPU:
      std::cout << " CPU:" << _cldevice.getInfo<CL_DEVICE_NAME>();
      break;
    case CL_DEVICE_TYPE_GPU:
      std::cout << " GPU:" << _cldevice.getInfo<CL_DEVICE_NAME>();
      break;
    default:
      std::cout << " DEFAULT:" << _cldevice.getInfo<CL_DEVICE_NAME>();
    }
  
  std::cout << std::endl;

  //Make a command queue
  _clcmdq = cl::CommandQueue(_clcontext, _cldevice);
}

void CLGLWindow::CallBackDisplayFunc(void)
{ 
  //Prepare for the OpenCL ticks
  glFinish();//Finish with the GL buffers
  //Setup the timings
  _currFrameTime = glutGet(GLUT_ELAPSED_TIME);

  //Run every objects OpenCL stage
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->clTick(_clcmdq, _clcontext);

  //Flush the OpenCL queue, so GL can use the buffers
  _clcmdq.finish();
  
  //Prepare for the GL render
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //Camera Positioning
  CameraSetup();

  //Enter the render ticks for all objects
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->glRender();

  drawAxis();
  
  drawArrow(_cameraDirection + Vector(-1,0,-1), Vector(-1,0,-1));

  
  glutSwapBuffers();

  ++frameCounter; 

  if (FPSmode && (_currFrameTime - _FPStime > 1000))
    {
      float fps = frameCounter * 1000.0 / (_currFrameTime - _FPStime);
      
      std::stringstream newWinTitle;
      newWinTitle << getWindowTitle() << " FPS : " << fps;
      
      glutSetWindowTitle(newWinTitle.str().c_str());
      frameCounter = 0;
      _FPStime = _currFrameTime;
    }

  _lastFrameTime = _currFrameTime;
}

void CLGLWindow::drawAxis()
{
  GLdouble nearPlane = 0.1,
    axisScale = 0.05;
  
  //We're drawing an overlayed axis so disable depth testing
  glViewport(0,0,100,100);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity();
  gluPerspective(45.0f, 1, 0.1f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity();  
  
  //near plane is at 0.1, the axis are 0.25 long so
  
  glTranslatef (0, 0, -(nearPlane + axisScale));

  glColor4f (4.0/256,104.0/256.0,202.0/256.0, 0.7); // Color the axis box a transparent blue
  glBegin(GL_QUADS);		
  glVertex3f(-1,-1, 0);
  glVertex3f( 1,-1, 0);
  glVertex3f( 1, 1, 0);
  glVertex3f(-1, 1, 0);
  glEnd();

  glRotatef(_rotatey, 1.0, 0.0, 0.0);
  glRotatef(_rotatex, 0.0, 1.0, 0.0);
  //glRotatef (tip , 1,0,0);
  //glRotatef (turn, 0,1,0);
  glScalef (axisScale, axisScale, axisScale);

  glLineWidth (2.0);
    
  glColor3f (1,0,0); // X axis is red.
  drawArrow(Vector(1,0,0), Vector(0,0,0));
  glColor3f (0,1,0); // Y axis is green.
  drawArrow(Vector(0,1,0), Vector(0,0,0));
  glColor3f (0,0,1); // Z axis is blue.
  drawArrow(Vector(0,0,1), Vector(0,0,0));
  
  //Do the text
  glColor3f(1,1,1);
  GLScribe::cout << GLScribe::cursor(1,0,0) << "X"
		 << GLScribe::cursor(0,1,0) << "Y"
		 << GLScribe::cursor(0,0,1) << "Z";

  glMatrixMode(GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix ();

  glViewport(0, 0, _width,_height);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
}

void CLGLWindow::CallBackReshapeFunc(int w, int h)
{
  _width = w;
  _height = h;
  
  //Setup the viewport
  glViewport(0, 0, _width, _height); 
  //Now reset the projection matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();


  //gluPerspective(45.0f, , 0.1f, 1000.0f);

  GLdouble fovy = 45.0f,
    zNear = 0.1f,
    zFar = 1000.0f,
    aspect = ((GLdouble)_width) / _height;
  GLdouble xmin, xmax, ymin, ymax;
  
  ymax = zNear * std::tan(fovy * M_PI / 360.0);
  ymin = -ymax;
  xmin = ymin * aspect;
  xmax = ymax * aspect;
  glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);

  glMatrixMode(GL_MODELVIEW);
}

void 
CLGLWindow::CallBackIdleFunc(void)
{
  CallBackDisplayFunc();
}

void 
CLGLWindow::setWindowtitle(const std::string& newtitle) 
{ 
  windowTitle = newtitle;
  glutSetWindowTitle(windowTitle.c_str());
}

void 
CLGLWindow::displayFPS(bool enable) 
{ 
  if (enable && !FPSmode) 
    {
      _FPStime = _currFrameTime; 
      frameCounter = 0;

      glutSetWindowTitle((windowTitle + " FPS : N/A").c_str());
      FPSmode = true;
    } 
  else if (!enable && FPSmode)
    {
      glutSetWindowTitle(windowTitle.c_str());
      FPSmode = false;
    }
}

void 
CLGLWindow::CallBackMouseFunc(int button, int state, int x, int y)
{
  switch (button)
    {
    case GLUT_LEFT_BUTTON:

      if (state == GLUT_DOWN)
	{
	  _oldMouseX = x;
	  _oldMouseY = y;
	  
	  keyState |= LEFTMOUSE;
	}
      else
	keyState &= ~LEFTMOUSE;
      break;
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN)
	{
	  _oldMouseX = x;
	  _oldMouseY = y;
	  
	  keyState |= RIGHTMOUSE;
	}
      else
	keyState &= ~RIGHTMOUSE;
      break;
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN)
	{
	  _oldMouseX = x;
	  _oldMouseY = y;
	  
	  keyState |= MIDDLEMOUSE;
	}
      else
	keyState &= ~MIDDLEMOUSE;
      break;
    default:
      break;
    }
}

void 
CLGLWindow::CallBackMouseWheelFunc(int button, int dir, int x, int y)
{
  //viewscale *= (dir > 0) ? 0.9 : (1/0.9);
}

void 
CLGLWindow::CallBackMotionFunc(int x, int y)
{
  float diffY = (y-_oldMouseY) * _mouseSensitivity;
  float diffX = (x-_oldMouseX) * _mouseSensitivity;

  switch (keyState)
    {
    case LEFTMOUSE:
      _rotatex += diffX;
      _rotatey = clamp(diffY + _rotatey, -90, 90);
      break;
    case RIGHTMOUSE:
      _cameraZ += (y-_oldMouseY) * _mouseSensitivity * 0.05;
      break;
    case MIDDLEMOUSE:
      _cameraX += (y-_oldMouseY) * _mouseSensitivity * 0.05;
      _cameraY += (x-_oldMouseX) * _mouseSensitivity * 0.05;
      break;
    default:
      break;
    }

  _oldMouseX = x;
  _oldMouseY = y;
}

void 
CLGLWindow::CallBackKeyboardFunc(unsigned char key, int x, int y)
{
  keyStates[key] = true;

  switch (key)
    {
      ///SPECIAL KEYPRESSES
    case 'F':
      displayFPS(false);
      break;
    case 'f':
      displayFPS(true);
    break;
    case 't':
    case 'T':
      for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
	   iPtr != RenderObjects.end(); ++iPtr)
	(*iPtr)->setRenderMode(RenderObj::TRIANGLES);
    break;
    case 'l':
    case 'L':
      for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
	   iPtr != RenderObjects.end(); ++iPtr)
	(*iPtr)->setRenderMode(RenderObj::LINES);
    break;
    case 'p':
    case 'P':
      for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
	   iPtr != RenderObjects.end(); ++iPtr)
	(*iPtr)->setRenderMode(RenderObj::POINTS);
    break;
    default:
      break;
    }
}

void 
CLGLWindow::CallBackKeyboardUpFunc(unsigned char key, int x, int y)
{
  keyStates[key] = false;
}

void
CLGLWindow::CallBackSpecialFunc(int key, int x, int y)
{
  //specialKeys = glutGetModifiers();
}

void
CLGLWindow::CallBackSpecialUpFunc(int key, int x, int y)
{
  //specialKeys = glutGetModifiers();  
}


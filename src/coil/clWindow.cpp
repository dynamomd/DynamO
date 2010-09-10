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

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <magnet/GLBuffer.hpp>
#include "glprimatives/glscribe.hpp"
#include "glprimatives/arrow.hpp"

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
  _mouseSensitivity(0.3),
  _moveSensitivity(0.00005),
  _specialKeys(0),
  _hostTransfers(hostTransfers),
  _shadows(false),
  _shadowMapSize(1024)
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
  _viewPortInfo._cameraZ -= _forward * moveAmp * std::cos(_viewPortInfo._rotatey * (M_PI/ 180)) 
    * std::sin(_viewPortInfo._rotatex  * (M_PI/ 180) + M_PI * 0.5);  
  _viewPortInfo._cameraX -= _forward * moveAmp * std::cos(_viewPortInfo._rotatey * (M_PI/ 180)) 
    * std::cos(_viewPortInfo._rotatex  * (M_PI/ 180) + M_PI * 0.5);
  _viewPortInfo._cameraY += -_forward * moveAmp * std::sin(_viewPortInfo._rotatey * (M_PI/ 180));

  //Strafe movement
  _viewPortInfo._cameraZ += _sideways * moveAmp * std::sin(_viewPortInfo._rotatex * (M_PI/ 180));
  _viewPortInfo._cameraX += _sideways * moveAmp * std::cos(_viewPortInfo._rotatex * (M_PI/ 180));

  //Vertical movement
  _viewPortInfo._cameraY += _vertical * moveAmp;

  glLoadIdentity();
  //gluLookAt(-viewscale, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  glRotatef(_viewPortInfo._rotatey, 1.0, 0.0, 0.0);
  glRotatef(_viewPortInfo._rotatex, 0.0, 1.0, 0.0);
  glTranslatef(-_viewPortInfo._cameraX,-_viewPortInfo._cameraY,-_viewPortInfo._cameraZ);
  
  //store the matricies for shadow calculations
  glGetFloatv(GL_MODELVIEW_MATRIX, _viewPortInfo._viewMatrix);
  glGetFloatv(GL_PROJECTION_MATRIX, _viewPortInfo._projectionMatrix);

  Matrix viewTransform = Rodrigues(Vector(0,-_viewPortInfo._rotatex * M_PI/180,0)) 
    * Rodrigues(Vector(-_viewPortInfo._rotatey * M_PI/180.0,0,0));

  _viewPortInfo._cameraDirection =  viewTransform * Vector(0,0,-1);
  _viewPortInfo._cameraUp = viewTransform * Vector(0,1,0);
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

  //Check for shadow support
  _shadows = true;
  if (!glewIsSupported("GL_ARB_depth_texture"))
    {
      std::cout << "GL_ARB_depth_texture not supported, shadows disabled\n";
      _shadows = false;
    }

  if (!glewIsSupported("GL_ARB_shadow"))
    {
      std::cout << "GL_ARB_shadow not supported, shadows disabled\n";
      _shadows = false;
    }
   
  //Now begins the example
  glClearColor(0.8f, 0.8f, 0.8f, 1.0f);

  glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);

  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  //Make OpenGL renormalize lighting vectors for us (incase we use glScale)
  glEnable(GL_NORMALIZE);

  //We need to cull for shadows
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW); //The default

  //Both the front and back materials track the current color
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL); //and enable it

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Blend colors using the alpha channel

  glShadeModel(GL_SMOOTH);

  //Shadow map initialisation
  if (_shadows)
    {
      glGenTextures(1, &_shadowMapTexture);
      glBindTexture(GL_TEXTURE_2D, _shadowMapTexture);
      glTexImage2D( GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, _shadowMapSize, _shadowMapSize, 0,
		    GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    }

  //Setup the viewport
  CallBackReshapeFunc(_width, _height);
  CameraSetup();

  //Light our scene!
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);   

  //Light number one
  //Position is set in the CameraSetup!
  //Ambient lighting
  GLfloat ambient_light[] = {0.0f, 0.0f, 0.0f, 1.0f}; 
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);

  _light0 = lightInfo(Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 0.0f), GL_LIGHT0);
  
  GLfloat specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  GLfloat specShininess[] = { 0.0f };
  GLfloat specular[] = {0.0, 0.0, 0.0, 1.0};
  glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
  glMaterialfv(GL_FRONT, GL_SHININESS, specShininess);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

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
  _clcmdq = cl::CommandQueue(_clcontext, _cldevice/*, CL_QUEUE_PROFILING_ENABLE*/) ;
}

void CLGLWindow::CallBackDisplayFunc(void)
{ 
  //Prepare for the OpenCL ticks
  glFinish();//Finish with the GL buffers
  //Setup the timings
  _currFrameTime = glutGet(GLUT_ELAPSED_TIME);

  const float rotateSpeed = 1000;
  _light0 = lightInfo(Vector(2.0f * std::cos(_currFrameTime/rotateSpeed), 
			     2.0f, 
			     2.0f * std::sin(_currFrameTime/rotateSpeed)), 
		      Vector(0.0f, 0.5f, 0.0f), GL_LIGHT0,
		      45.0f, 5, 0.01f);


  //Run every objects OpenCL stage
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->clTick(_clcmdq, _clcontext);

  //Camera Positioning
  CameraSetup();

  //Flush the OpenCL queue, so GL can use the buffers
  _clcmdq.finish();
  
  //Prepare for the GL render
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //If shadows are enabled, we must draw first from the lights perspective
  if (_shadows)
    {
      VECTOR4D white(1,1,1,1);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, white*0.3f);
      glLightfv(GL_LIGHT0, GL_AMBIENT, white*0.3f);

      ////Store the camera matrices and load the light's
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadMatrixf(_light0._projectionMatrix);
      
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadMatrixf(_light0._viewMatrix);

      //Now render the scene from the lights perspective
      //The viewport should change to the shadow maps size
      glViewport(0, 0, _shadowMapSize, _shadowMapSize);

      //Draw back faces into the shadow map
      glCullFace(GL_FRONT);
      
      //Disable color writes, and use flat shading for speed
      glShadeModel(GL_FLAT);
      glColorMask(0, 0, 0, 0);
            
      //Render the scene
      drawScene();
      
      //Copy the shadow texture
      glBindTexture(GL_TEXTURE_2D, _shadowMapTexture);
      glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, _shadowMapSize, _shadowMapSize);

      //Restore the draw mode
      glCullFace(GL_BACK);
      glShadeModel(GL_SMOOTH);
      glColorMask(1, 1, 1, 1);
      
      //Restore the viewport
      glViewport(0, 0, _width,_height);
      
      ////Restore the Camera matricies
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();

      //////////////////Pass 2//////////////////
      //Only clear the depth buffer, the color buffer is already clear
      glClear(GL_DEPTH_BUFFER_BIT); 

      drawScene();

      //////////////////Pass 3//////////////////
      //Setup a bright light
      glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
      glLightfv(GL_LIGHT0, GL_SPECULAR, white);

      //Now setup the texgen matrix
      static MATRIX4X4 biasMatrix(0.5f, 0.0f, 0.0f, 0.0f,
				  0.0f, 0.5f, 0.0f, 0.0f,
				  0.0f, 0.0f, 0.5f, 0.0f,
				  0.5f, 0.5f, 0.5f, 1.0f); //bias from [-1, 1] to [0, 1]

      MATRIX4X4 textureMatrix = biasMatrix 
	* _light0._projectionMatrix * _light0._viewMatrix;
      
      //Set up texture coordinate generation.
      glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
      glTexGenfv(GL_S, GL_EYE_PLANE, textureMatrix.GetRow(0));
      
      glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
      glTexGenfv(GL_T, GL_EYE_PLANE, textureMatrix.GetRow(1));
      
      glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
      glTexGenfv(GL_R, GL_EYE_PLANE, textureMatrix.GetRow(2));
      
      glTexGeni(GL_Q, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
      glTexGenfv(GL_Q, GL_EYE_PLANE, textureMatrix.GetRow(3));

      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
      glEnable(GL_TEXTURE_GEN_R);
      glEnable(GL_TEXTURE_GEN_Q);

      //Bind & enable shadow map texture
      glBindTexture(GL_TEXTURE_2D, _shadowMapTexture);
      glEnable(GL_TEXTURE_2D);
      
      //Enable shadow comparison
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE);
      
      //Shadow comparison should be true (ie not in shadow) if r<=texture
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);
      
      //Shadow comparison should generate an INTENSITY result
      glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE_ARB, GL_INTENSITY);
      
      //Set alpha test to discard false comparisons
      glAlphaFunc(GL_GEQUAL, 0.99f);
      glEnable(GL_ALPHA_TEST);

      drawScene();

      //Disable textures and texgen
      glDisable(GL_TEXTURE_2D);
      
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_GEN_R);
      glDisable(GL_TEXTURE_GEN_Q);
      
      //Restore other states
      glDisable(GL_ALPHA_TEST);
    }
  else
    drawScene();

  drawAxis();

  //Draw the light source
  glColor3f(1,1,0);
  glPushMatrix();
  glTranslatef(_light0._position.x, _light0._position.y, _light0._position.z);
  glutSolidSphere(0.1f, 5, 5);
  glPopMatrix();


  //coil::glprimatives::drawArrow(_cameraDirection + Vector(-1,0,-1), Vector(-1,0,-1)); 
  //coil::glprimatives::drawArrow(_cameraUp + Vector(-1,0,-1), Vector(-1,0,-1));
 
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

void 
CLGLWindow::drawScene()
{
  //Move the world lights
  GLfloat light0_position[] = {_light0._position.x, _light0._position.y, _light0._position.z, 0.0f};
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  //Enter the render ticks for all objects
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->glRender();

  //Draw a ground
  glColor3f(1,1,1);
  glPushMatrix();
  glTranslatef(0.0f,-0.6f,0.0f);
  glScalef(4.0f, 0.01f, 4.0f);
  glutSolidCube(1.0);

  glPopMatrix();

}


void CLGLWindow::drawAxis()
{
  GLdouble nearPlane = 0.1,
    axisScale = 0.07;
  
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

  glEnable(GL_BLEND); //Enable blending

  glColor4f (4.0/256,104.0/256.0,202.0/256.0, 0.7); // Color the axis box a transparent blue
  glBegin(GL_QUADS);		
  glVertex3f(-1,-1, 0);
  glVertex3f( 1,-1, 0);
  glVertex3f( 1, 1, 0);
  glVertex3f(-1, 1, 0);
  glEnd();

  glDisable(GL_BLEND); //Turn blending back off

  glRotatef(_viewPortInfo._rotatey, 1.0, 0.0, 0.0);
  glRotatef(_viewPortInfo._rotatex, 0.0, 1.0, 0.0);
  glScalef (axisScale, axisScale, axisScale);

  glLineWidth (2.0);
    
  glColor3f (1,0,0); // X axis is red.
  coil::glprimatives::drawArrow(Vector( 0.5,-0.5,-0.5), 
				Vector(-0.5,-0.5,-0.5));

  glColor3f (0,1,0); // Y axis is green.
  coil::glprimatives::drawArrow(Vector(-0.5, 0.5,-0.5), 
				Vector(-0.5,-0.5,-0.5));

  glColor3f (0,0,1); // Z axis is blue.
  coil::glprimatives::drawArrow(Vector(-0.5,-0.5, 0.5), 
				Vector(-0.5,-0.5,-0.5));
  
  //Do the text
  glColor3f(1,1,1);
  GLScribe::cout << GLScribe::cursor( 0.5,-0.5,-0.5) << "X"
		 << GLScribe::cursor(-0.5, 0.5,-0.5) << "Y"
		 << GLScribe::cursor(-0.5,-0.5, 0.5) << "Z";

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

  _viewPortInfo._aspectRatio = ((GLdouble)_width) / _height;

  gluPerspective(_viewPortInfo._fovY, _viewPortInfo._aspectRatio, _viewPortInfo._zNearDist, 
		 _viewPortInfo._zFarDist);

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
  if (dir > 0)
    _moveSensitivity *= 1.1;
  else
    _moveSensitivity /= 1.1;

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
      _viewPortInfo._rotatex += diffX;
      _viewPortInfo._rotatey = clamp(diffY + _viewPortInfo._rotatey, -90, 90);
      break;
    case RIGHTMOUSE:
      _viewPortInfo._cameraZ += (y-_oldMouseY) * _mouseSensitivity * 0.05;
      break;
    case MIDDLEMOUSE:
      _viewPortInfo._cameraX += (y-_oldMouseY) * _mouseSensitivity * 0.05;
      _viewPortInfo._cameraY += (x-_oldMouseX) * _mouseSensitivity * 0.05;
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


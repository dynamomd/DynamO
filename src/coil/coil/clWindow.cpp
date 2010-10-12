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

#include <locale>

inline float clamp(float x, float a, float b)
{
  return x < a ? a : (x > b ? b : x);
}

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "glprimatives/glscribe.hpp"
#include "glprimatives/arrow.hpp"

#include <magnet/CL/CLGL.hpp>

CLGLWindow::CLGLWindow(int setWidth, int setHeight,
                       int initPosX, int initPosY,
                       std::string title
		       ):
  _height(setHeight),
  _width(setWidth),
  _windowX(initPosX),
  _windowY(initPosY),
  keyState(DEFAULT),
  windowTitle(title),
  _frameCounter(0),
  _updateCounter(0),
  _mouseSensitivity(0.3),
  _moveSensitivity(0.001),
  _specialKeys(0),
  _shaderPipeline(false)
{
  for (size_t i(0); i < 256; ++i) keyStates[i] = false;
}

CLGLWindow::~CLGLWindow()
{
  deinit(true);
}

void 
CLGLWindow::initOpenGL()
{
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowSize(_width, _height);
  glutInitWindowPosition(_windowX, _windowY);

  CoilMaster::getInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

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
  _shaderPipeline = true;
  if (!GLEW_ARB_depth_texture || !GLEW_ARB_shadow)
    {
      std::cout << "GL_ARB_depth_texture or GL_ARB_shadow not supported.\n";
      _shaderPipeline = false;
    }
  else if (!GLEW_ARB_fragment_program || !GLEW_ARB_vertex_program
	   || !GLEW_ARB_fragment_shader || !GLEW_ARB_vertex_shader)
    {
      std::cout << "OpenGL driver doesn't support programmable shaders.\n";
      _shaderPipeline = false;
    }

  if (!_shaderPipeline)
    std::cout << "Shader pipeline disabled.\n"
	      << "This also disables all other effects.\n";

  glDrawBuffer(GL_BACK);

  glClearColor(0.8f, 0.8f, 0.8f, 1.0f);

  glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);

  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  //We need blending
  glEnable(GL_BLEND);
  //Blend colors using the alpha channel
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

  //Make OpenGL renormalize lighting vectors for us (incase we use glScale)
  glEnable(GL_NORMALIZE);

  //Switch on line aliasing
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  //We need to cull for shadows
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW); //The default

  //Both the front and back materials track the current color
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL); //and enable it


  glShadeModel(GL_SMOOTH);

  //Setup the viewport
  CallBackReshapeFunc(_width, _height);
  _viewPortInfo.CameraSetup();

  //Light our scene!
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);   

  //Light number one
  //Position is set in the CameraSetup!
  //Ambient lighting
  GLfloat ambient_light[] = {0.0f, 0.0f, 0.0f, 1.0f}; 
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);


  _light0 = magnet::GL::lightInfo(Vector(1.5f,  1.5f, 1.0f),//Position
				  Vector(0.0f, -0.3f, 0.0f),//Lookat
				  GL_LIGHT0, //GL handle
				  45.0f,//Beam angle
				  50,//rangeMax
				  0.005//rangeMin
				  );
  
  GLfloat specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  GLfloat specShininess[] = { 0.0f };
  GLfloat specular[] = {0.0, 0.0, 0.0, 1.0};
  glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
  glMaterialfv(GL_FRONT, GL_SHININESS, specShininess);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

  //Setup the keyboard controls
  glutIgnoreKeyRepeat(1);

  _FPStime = glutGet(GLUT_ELAPSED_TIME);

  //Build the offscreen rendering FBO's
  if (_shaderPipeline)
    {
      _shadowFBO.init(1024);
      _shadowShader.build();
    }

  //Now init the render objects  
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initOpenGL();

  if (_shaderPipeline)
    {
      _renderTarget.reset(new magnet::GL::FBO());
      _renderTarget->init(_width, _height);
    }
}

void 
CLGLWindow::initOpenCL()
{
  _CLState.init();
  
  if (cl::GLBuffer::hostTransfers()) 
    std::cout << "\n!!!!!!!Host transfers have been enabled!!!!!!, slow performance is expected\n";

  //Now init the render objects  
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initOpenCL(_CLState);
}

//The glade xml file is "linked" into a binary file and stuffed in the executable, these are the symbols to its data
extern const char _binary_src_coil_coil_clwingtk_gladexml_start[];
extern const char _binary_src_coil_coil_clwingtk_gladexml_end[];

void
CLGLWindow::initGTK()
{
  {
    Glib::ustring glade_data
      (reinterpret_cast<const char *>(_binary_src_coil_coil_clwingtk_gladexml_start), 
       _binary_src_coil_coil_clwingtk_gladexml_end
       -_binary_src_coil_coil_clwingtk_gladexml_start);
    
    _refXml = Gtk::Builder::create_from_string(glade_data);
  }

  _timeout_connection
    = Glib::signal_timeout().connect_seconds(sigc::mem_fun(this, &CLGLWindow::GTKTick), 2);

  _refXml->get_widget("controlWindow", controlwindow);
  
  if (_shaderPipeline)
    {
      {//Enable the whole shader frame
	Gtk::Frame* shaderFrame;
	_refXml->get_widget("RenderPipelineFrame", shaderFrame);
	shaderFrame->set_sensitive(true);
      }
      
      {//Setup the checkbox
	Gtk::CheckButton* shaderEnable;
	_refXml->get_widget("ShaderPipelineEnable", shaderEnable);
	
	shaderEnable->set_active(true);

	shaderEnable->signal_toggled().connect(sigc::mem_fun(this, &CLGLWindow::pipelineEnableCallback));
      }

      if (GLEW_EXT_framebuffer_multisample)
	{//Offer anti aliasing
	  {//Turn on the antialiasing box
	    Gtk::HBox* multisampleBox;
	    _refXml->get_widget("multisampleBox", multisampleBox);
	    multisampleBox->set_sensitive(true);
	  }

	  {//Connect the anti aliasing checkbox
	    Gtk::CheckButton* multisampleEnable;
	    _refXml->get_widget("multisampleEnable", multisampleEnable);
	    multisampleEnable->signal_toggled()
	      .connect(sigc::mem_fun(*this, &CLGLWindow::multisampleEnableCallback));
	  }
	  
	  Gtk::ComboBox* aliasSelections;
	  _refXml->get_widget("multisampleLevels", aliasSelections);

	  struct aliasColumns : public Gtk::TreeModel::ColumnRecord
	  {
	    aliasColumns() { add(m_col_id); }
	    Gtk::TreeModelColumn<int> m_col_id;
	  };

	  aliasColumns vals;
	  Glib::RefPtr<Gtk::ListStore> m_refTreeModel = Gtk::ListStore::create(vals);
	  aliasSelections->set_model(m_refTreeModel);

	  GLint maxSamples;
	  glGetIntegerv(GL_MAX_SAMPLES, &maxSamples);

	  for ( ; maxSamples > 1; maxSamples >>= 1)
	    {
	      Gtk::TreeModel::Row row = *(m_refTreeModel->prepend());
	      row[vals.m_col_id] = maxSamples;
	    }

	  aliasSelections->pack_start(vals.m_col_id);
	  aliasSelections->set_active(0);

	  aliasSelections->signal_changed()
	    .connect(sigc::mem_fun(*this, &CLGLWindow::multisampleEnableCallback));
	}
    }
}

bool 
CLGLWindow::GTKTick()
{
  //This callback is used to calculate the FPS and sim update rates
  int currFrameTime = glutGet(GLUT_ELAPSED_TIME);

  float fps = _frameCounter * 1000.0 / (currFrameTime - _FPStime);
  float ups = _updateCounter * 1000.0 / (currFrameTime - _FPStime);
  
  std::stringstream fpsstring;
  fpsstring << "FPS:" << fps;
  
  Gtk::Label* label;
  _refXml->get_widget("RenderUpdateLabel", label);
  label->set_text(fpsstring.str());

  fpsstring.str("");
  fpsstring << "UPS:" << ups;

  _refXml->get_widget("SimUpdateLabel", label);
  label->set_text(fpsstring.str());

  _frameCounter = 0;
  _updateCounter = 0;
  _FPStime = currFrameTime;

  return true;
}

void 
CLGLWindow::pipelineEnableCallback()
{
  Gtk::CheckButton* shaderEnable;
  _refXml->get_widget("ShaderPipelineEnable", shaderEnable);
  
  _shaderPipeline = shaderEnable->get_active();
}

void 
CLGLWindow::multisampleEnableCallback()
{
  Gtk::CheckButton* multisampleEnable;
  _refXml->get_widget("multisampleEnable", multisampleEnable);
  if (multisampleEnable->get_active())
    {
      Gtk::ComboBox* aliasSelections;
      _refXml->get_widget("multisampleLevels", aliasSelections);

      _renderTarget.reset(new magnet::GL::multisampledFBO(2 << aliasSelections->get_active_row_number()));
      _renderTarget->init(_width, _height);
    }
  else
    {
      _renderTarget.reset(new magnet::GL::FBO());
      _renderTarget->init(_width, _height);
    }
}

void
CLGLWindow::init()
{
  magnet::thread::ScopedLock lock(_destroyLock);

  if (_readyFlag) return;

  initOpenGL();
  initOpenCL();
  initGTK();
  _readyFlag = true;
}


void
CLGLWindow::deinit(bool andGlutDestroy)
{
  magnet::thread::ScopedLock lock(_destroyLock);
  
  if (!_readyFlag) return;
  _readyFlag = false;

  ////////////////////GTK

  _timeout_connection.disconnect();

  {
    Gtk::Window* controlwindow;
    _refXml->get_widget("controlWindow", controlwindow);  
    controlwindow->hide_all();
  }
  
  _refXml.reset(); //Destroy GTK instance

  /////////////////OpenCL
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    delete *iPtr;

  _CLState = magnet::CL::CLGLState();

  ///////////////////OpenGL
  if (_shaderPipeline)
    {
      _shadowFBO = magnet::GL::shadowFBO();
      _shadowShader = magnet::GL::shadowShader();
    }

  _renderTarget.reset();

  ///////////////////Finally, unregister with COIL
  CoilMaster::getInstance().CallGlutDestroyWindow(this, andGlutDestroy);
}

void 
CLGLWindow::CallBackDisplayFunc(void)
{
  if (!CoilMaster::getInstance().isRunning()) return;

  //Prepare for the OpenCL ticks
  glFinish();//Finish with the GL buffers
  //Setup the timings
  int _currFrameTime = glutGet(GLUT_ELAPSED_TIME);

//  const float speed = 1000;
//  _light0 = lightInfo(Vector(1.5f*std::cos(_currFrameTime/speed), 1.5f, 
//			     1.5f * std::sin(_currFrameTime/speed)), 
//		      Vector(0.0f, 0.0f, 0.0f), GL_LIGHT0);

  //Run every objects OpenCL stage
  for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->clTick(_CLState);

  //Camera Positioning

  float moveAmp = (_currFrameTime - _lastFrameTime) * _moveSensitivity;      
  float forward  = moveAmp * ( keyStates['w'] - keyStates['s']);
  float sideways = moveAmp * ( keyStates['d'] - keyStates['a']);
  float vertical = moveAmp * ( keyStates['q'] - keyStates['z']);
  _viewPortInfo.CameraSetup(forward, sideways, vertical);

  //Flush the OpenCL queue, so GL can use the buffers
  _CLState.getCommandQueue().finish();
  
  //Prepare for the GL render
  if (_shaderPipeline)
    {
      //////////////////Pass 1//////////////////
      ///Here we draw from the 
      _shadowFBO.setup(_light0);

#ifdef GL_VERSION_1_1
      glEnable (GL_POLYGON_OFFSET_FILL);
      glPolygonOffset (1., 1.);
#endif 
      
      drawScene();
      
#ifdef GL_VERSION_1_1
      glDisable (GL_POLYGON_OFFSET_FILL);
#endif

      _shadowFBO.restore();
      //////////////////Pass 2//////////////////
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

      //In both cases we use the texture matrix, instead of the EYE_PLANE
      //We bind to the 7th texture unit???? 
      glActiveTextureARB(GL_TEXTURE7);
      glMatrixMode(GL_TEXTURE);

      _light0.buildShadowTextureMatrix();

      MATRIX4X4 invView = _viewPortInfo._viewMatrix.GetInverse();
      glMultMatrixf(invView);
	
      glMatrixMode(GL_MODELVIEW);

      //Bind to the multisample buffer
      _renderTarget->attach();
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

      //Setup and run the shadow shader
      glBindTexture(GL_TEXTURE_2D, _shadowFBO.getShadowTexture());
      _shadowShader.attach(_shadowFBO.getShadowTexture(), _shadowFBO.getLength(), 7);
      drawScene();
      
      _renderTarget->detach();

      //Restore the fixed pipeline
      //And turn off the shadow texture
      glUseProgramObjectARB(0);

      //Now blit the stored scene to the screen
      _renderTarget->blitToScreen(_width, _height);

      /////////////FILTERING
      //Now we blur the output from the offscreen FBO
      //
      //glActiveTextureARB(GL_TEXTURE0);
      //
      //_FBO1.attach();
      //glBindTexture(GL_TEXTURE_2D, _myFBO.getColorTexture());
      ////_laplacianFilter.invoke(0, _width, _height);
      //_blurFilter.invoke(0, _width, _height);
      //
      //_FBO1.detach();
      //_FBO1.blitToScreen(_width, _height);      
    }
  else    
    {      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);       
      drawScene();
    }
  
  glClear(GL_DEPTH_BUFFER_BIT);       
  drawAxis();

  //Draw the light source
  _light0.drawLight();

  glutSwapBuffers();

  ++_frameCounter; 
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
  
   glBegin(GL_QUADS);
   //Front
   glNormal3f(0, 1, 0);
   glVertex3f(-100 + _viewPortInfo._cameraX, -0.51, -100 + _viewPortInfo._cameraZ);
   glVertex3f(-100 + _viewPortInfo._cameraX, -0.51,  100 + _viewPortInfo._cameraZ);
   glVertex3f( 100 + _viewPortInfo._cameraX, -0.51,  100 + _viewPortInfo._cameraZ);
   glVertex3f( 100 + _viewPortInfo._cameraX, -0.51, -100 + _viewPortInfo._cameraZ);

//   glNormal3f(0, -1, 0);
//   glVertex3f(-10 + _viewPortInfo._cameraX, - 0.52, -10 + _viewPortInfo._cameraZ);
//   glVertex3f(-10 + _viewPortInfo._cameraX, - 0.52,  10 + _viewPortInfo._cameraZ);
//   glVertex3f( 10 + _viewPortInfo._cameraX, - 0.52,  10 + _viewPortInfo._cameraZ);
//   glVertex3f( 10 + _viewPortInfo._cameraX, - 0.52, -10 + _viewPortInfo._cameraZ);

   //Back
//   glNormal3f(0, -1, 0);
//   glVertex3f(-10, -0.51, -10);
//   glVertex3f(-10, -0.51,  10);
//   glVertex3f( 10, -0.51,  10);
//   glVertex3f( 10, -0.51, -10);

   glEnd();
  
//  glPushMatrix();
//  glTranslatef(0.0f,-0.6f,0.0f);
//  glScalef(4.0f, 0.01f, 4.0f);
//  glutSolidCube(1.0);

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

  glColor4f (4.0/256,104.0/256.0,202.0/256.0, 0.7); // Color the axis box a transparent blue
  glBegin(GL_QUADS);		
  glVertex3f(-1,-1, 0);
  glVertex3f( 1,-1, 0);
  glVertex3f( 1, 1, 0);
  glVertex3f(-1, 1, 0);
  glEnd();

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
  if (!CoilMaster::getInstance().isRunning() || !_readyFlag) return;

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

  if (_shaderPipeline)
    _renderTarget->resize(_width, _height);
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
  keyStates[std::tolower(key)] = true;

  switch (key)
    {
      ///SPECIAL KEYPRESSES
    case 't':
      for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
	   iPtr != RenderObjects.end(); ++iPtr)
	(*iPtr)->setRenderMode(RenderObj::TRIANGLES);
      break;
    case 'l':
      for (std::vector<RenderObj*>::iterator iPtr = RenderObjects.begin();
	   iPtr != RenderObjects.end(); ++iPtr)
	(*iPtr)->setRenderMode(RenderObj::LINES);
      break;
    case 'p':
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
  keyStates[std::tolower(key)] = false;
}


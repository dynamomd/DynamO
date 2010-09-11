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

#include <vector>

#include <coil/glutMaster.hpp>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#include <CL/cl.hpp>

#include <magnet/GLBuffer.hpp>

#include <coil/extcode/static_assert.hpp>
#include <coil/extcode/vector2.hpp>
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/Maths/Maths.h>

#include <limits>

class CLGLWindow : public CoilWindow
{
public:
  CLGLWindow(int setWidth, int setHeight,
	     int setInitPositionX, int setInitPositionY,
	     std::string title,
	     cl::Platform& plat,
	     bool hostTransfers = false);

  ~CLGLWindow();
  
  virtual void CallBackDisplayFunc(void); 
  virtual void CallBackIdleFunc(void);
  void CallBackReshapeFunc(int w, int h);    

  cl::Platform& getCLPlatform() { return _clplatform; }
  cl::Context& getCLContext() { return  _clcontext; }
  cl::Device& getCLDevice() { return  _cldevice; }
  cl::CommandQueue& getCLCommandQueue() { return  _clcmdq; }

  const std::string& getWindowTitle() const { return windowTitle; }
  void setWindowtitle(const std::string& newtitle);
  
  void displayFPS(bool enable = true);

  void addRenderObj(RenderObj* nObj) { RenderObjects.push_back(nObj); }

  bool HostTransferModeAllowed() { return _hostTransfers; }
  
  template<class T>  
  void addRenderObj() {STATIC_ASSERT(false,'Invalid Render Object'); }
  template<class T, class T1> 
  void addRenderObj(T1) { STATIC_ASSERT(false,'Check Arg Types');}
  template<class T, class T1, class T2> 
  void addRenderObj(T1, T2) { STATIC_ASSERT(false,'Check Arg Types'); }
  template<class T, class T1, class T2, class T3> 
  void addRenderObj(T1, T2, T3) { STATIC_ASSERT(false,'Check Arg Types'); } 
  template<class T, class T1, class T2, class T3, class T4> 
  void addRenderObj(T1, T2, T3, T4) { STATIC_ASSERT(false,'Check Arg Types'); } 
  template<class T, class T1, class T2, class T3, class T4, class T5> 
  void addRenderObj(T1, T2, T3, T4, T5) { STATIC_ASSERT(false,'Check Arg Types'); } 
  template<class T, class T1, class T2, class T3, class T4, class T5, class T6> 
  void addRenderObj(T1, T2, T3, T4, T5, T6) { STATIC_ASSERT(false,'Check Arg Types'); } 

  struct viewPortInfoType
  {
    viewPortInfoType():
      _rotatex(180),
      _rotatey(0),
      _cameraX(0),
      _cameraY(0),
      _cameraZ(0),
      _fovY(45),
      _aspectRatio(1),
      _zNearDist(0.0001),
      _zFarDist(10)
    {}

    float _rotatex;
    float _rotatey;
    float _cameraX;
    float _cameraY;
    float _cameraZ;
    
    GLdouble _fovY;
    GLdouble _aspectRatio;
    GLdouble _zNearDist;
    GLdouble _zFarDist;
    
    Vector _cameraDirection, _cameraUp;

    MATRIX4X4 _projectionMatrix;
    MATRIX4X4 _viewMatrix;
  };
 
  struct lightInfo
  {
    lightInfo() {}

    lightInfo(Vector position, Vector lookAtPoint, GLenum lightHandle, 
	      GLfloat beamAngle = 45.0f,
	      GLfloat rangeMax = 2.0f, GLfloat rangeMin = 0.001):
      _position(position),
      _lookAtPoint(lookAtPoint),
      _lightHandle(lightHandle)
    {
      //Build the view matrix and so on
      glPushMatrix();

      glLoadIdentity();
      gluPerspective(beamAngle, 1.0f, rangeMin, rangeMax);
      glGetFloatv(GL_MODELVIEW_MATRIX, _projectionMatrix);

      glLoadIdentity();
      Vector directionNorm = (lookAtPoint - position);
      directionNorm /= directionNorm.nrm();

      GLfloat rotationAngle = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | directionNorm);


      Vector RotationAxis = Vector(0,0,-1) ^ directionNorm;
      float norm = RotationAxis.nrm();
      RotationAxis /= norm;
      if (norm < std::numeric_limits<double>::epsilon())
	RotationAxis = Vector(1,0,0);

      glRotatef(-rotationAngle, RotationAxis.x,RotationAxis.y,RotationAxis.z);
      glTranslatef(-position.x,-position.y,-position.z);
      glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);

      glPopMatrix();
    }

    Vector _position;
    Vector _lookAtPoint;

    MATRIX4X4 _projectionMatrix;
    MATRIX4X4 _viewMatrix;

    GLenum _lightHandle;
  };

protected:
  cl::Platform _clplatform;
  cl::Context _clcontext;
  cl::Device _cldevice;
  cl::CommandQueue _clcmdq;

  size_t _height, _width;

  std::vector<RenderObj*> RenderObjects;

  void CallBackSpecialUpFunc(int key, int x, int y);
  void CallBackSpecialFunc(int key, int x, int y);
  void CallBackKeyboardFunc(unsigned char key, int x, int y);
  void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
  void CallBackMouseWheelFunc(int button, int dir, int x, int y);
  void CallBackMouseFunc(int button, int state, int x, int y);
  void CallBackMotionFunc(int x, int y);

private:
  void CameraSetup();

  void initOpenGL(int initPosX, int initPosY);
  void initOpenCL();

  void drawAxis();

  void drawScene();

  enum KeyStateType
    {
      DEFAULT = 0,
      LEFTMOUSE = 1,
      RIGHTMOUSE = 2,
      MIDDLEMOUSE = 4
    };

  int keyState;
  
  std::string windowTitle;
  bool FPSmode;
  size_t frameCounter;


  int _currFrameTime;
  int _lastFrameTime;
  int _FPStime; 

  viewPortInfoType _viewPortInfo;
  
  
  bool keyStates[256];

  float _mouseSensitivity; 
  float _moveSensitivity;
 
  int _oldMouseX, _oldMouseY;
  int _specialKeys;

  bool _hostTransfers;

  bool _shadows;
  GLuint _shadowMapTexture;
  int _shadowMapSize;

  lightInfo _light0;
};

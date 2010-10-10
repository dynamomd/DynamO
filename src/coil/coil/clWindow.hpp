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

#include <gtkmm.h>

#include <coil/coilMaster.hpp>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#include <CL/cl.hpp>

#include <coil/Maths/Maths.h>
#include <coil/extcode/vector2.hpp>
#include <magnet/static_assert.hpp>

#include <magnet/GL/light.hpp>
#include <magnet/CL/GLBuffer.hpp>
#include <magnet/GL/shadowShader.hpp>
#include <magnet/GL/shadowFBO.hpp>
#include <magnet/GL/viewPort.hpp>
#include <magnet/GL/multisampledFBO.hpp>

#include <magnet/GL/blur.hpp>
#include <magnet/GL/laplacianFilter.hpp>

#include <coil/RenderObj/RenderObj.hpp>


class CLGLWindow : public CoilWindow
{
public:
  CLGLWindow(int setWidth, int setHeight,
	     int setInitPositionX, int setInitPositionY,
	     std::string title);

  ~CLGLWindow();
  
  virtual void CallBackDisplayFunc(void); 
  virtual void CallBackIdleFunc(void);
  void CallBackReshapeFunc(int w, int h);    

  const std::string& getWindowTitle() const { return windowTitle; }
  void setWindowtitle(const std::string& newtitle);
  
  void displayFPS(bool enable = true);

  void addRenderObj(RenderObj* nObj) { RenderObjects.push_back(nObj); }

  template<class T>  
  T& addRenderObj() {M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  template<class T, class T1> 
  T& addRenderObj(T1) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types);}

  template<class T, class T1, class T2> 
  T& addRenderObj(T1, T2) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  template<class T, class T1, class T2, class T3> 
  T& addRenderObj(T1, T2, T3) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  template<class T, class T1, class T2, class T3, class T4> 
  T& addRenderObj(T1, T2, T3, T4) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  template<class T, class T1, class T2, class T3, class T4, class T5> 
  T& addRenderObj(T1, T2, T3, T4, T5) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  template<class T, class T1, class T2, class T3, class T4, class T5, class T6> 
  T& addRenderObj(T1, T2, T3, T4, T5, T6) { M_STATIC_ASSERT(!sizeof(T),Check_Arg_Types); }

  inline const int& getLastFrameTime() const { return _lastFrameTime; }

  magnet::CL::CLGLState& getCLState() { return _CLState; }

  void init();
  void deinit();
  bool acquire();
  void release();

protected:
  magnet::GL::shadowShader _shadowShader;
  magnet::GL::shadowFBO _shadowFBO;

  //Frame buffers to flip flop data between
  magnet::GL::multisampledFBO _FBO1;

  size_t _height, _width;
  int _windowX, _windowY;

  std::vector<RenderObj*> RenderObjects;

  void CallBackSpecialUpFunc(int key, int x, int y) {}
  void CallBackSpecialFunc(int key, int x, int y) {} 
  void CallBackKeyboardFunc(unsigned char key, int x, int y);
  void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
  void CallBackMouseWheelFunc(int button, int dir, int x, int y);
  void CallBackMouseFunc(int button, int state, int x, int y);
  void CallBackMotionFunc(int x, int y);

private:
  magnet::thread::Mutex _destroyLock;

  void CameraSetup();

  magnet::CL::CLGLState _CLState;

  
  virtual void initOpenGL();
  virtual void initOpenCL();

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

  magnet::GL::viewPort _viewPortInfo;
    
  bool keyStates[256];

  float _mouseSensitivity; 
  float _moveSensitivity;
 
  int _oldMouseX, _oldMouseY;
  int _specialKeys;

  bool _shaderPipeline;

  magnet::GL::lightInfo _light0;

  /////////GTK members
  virtual void initGTK();
  Glib::RefPtr<Gtk::Builder> _refXml;
};

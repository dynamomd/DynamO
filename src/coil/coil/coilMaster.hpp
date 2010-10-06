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

#include <gtkmm.h>

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <coil/coilWindow.hpp>
#include <magnet/thread/mutex.hpp>
#include <magnet/thread/taskQueue.hpp>

#include <map>
#include <vector>

class CoilMaster{
public:

  inline static 
  CoilMaster& getInstance()
  {
    static CoilMaster singletonInstance;
    return singletonInstance;
  }
  
  //Only for window objects to call
  void  CallGlutCreateWindow(const char*, CoilWindow*);
  void  CallGlutDestroyWindow(CoilWindow*);

  void addWindow(CoilWindow* window)
  {
    _renderQueue.queueTask(magnet::function::Task::makeTask(addWindowFunc, window));

    //Spinlock waiting for the window to initialize
    while (!window->isReady()) { smallSleep(); }
  }

  void bootCoil();

  inline void 
  shutdownCoil() 
  { 
    _runFlag = false; 
    _renderReadyFlag = false;
    _windowReadyFlag = false;
  }

  void waitForShutdown();

  inline bool isRunning() { return _runFlag; }

private:
  CoilMaster();
  ~CoilMaster();

  static inline void addWindowFunc(CoilWindow* window)
  { 
    window->init();
  }


  void smallSleep();

  volatile bool _runFlag; 
  volatile bool _renderReadyFlag;
  volatile bool _windowReadyFlag;
   
  ///////////////////////////Glut GL render layer//////////////////////////
  magnet::thread::Thread _renderThread;
  magnet::thread::TaskQueue _renderQueue;

  std::map<int, CoilWindow*> _viewPorts;
  std::map<CoilWindow*, int> _windows;
  static void CallBackDisplayFunc(); 
  static void CallBackCloseWindow();
  static void CallBackKeyboardFunc(unsigned char key, int x, int y);
  static void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
  static void CallBackMotionFunc(int x, int y);
  static void CallBackMouseFunc(int button, int state, int x, int y);
  static void CallBackMouseWheelFunc(int button, int dir, int x, int y);
  static void CallBackPassiveMotionFunc(int x, int y);
  static void CallBackReshapeFunc(int w, int h); 
  static void CallBackSpecialFunc(int key, int x, int y);   
  static void CallBackSpecialUpFunc(int key, int x, int y);   
  static void CallBackVisibilityFunc(int visible);
    
  void renderThreadEntryPoint();

  ///////////////////////////GTK window layer/////////////////////////////
  void windowThreadEntryPoint();
  bool onControlWindowDelete(GdkEventAny * pEvent);

  Gtk::Main _GTKit;
  Glib::RefPtr<Gtk::Builder> _refXml;

  bool GTKIldeFunc();

  magnet::thread::Thread _windowThread;
  magnet::thread::TaskQueue _windowQueue;
};

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
    if (!isRunning()) M_throw() << "Coil is not running, cannot add a window";

    _coilQueue.queueTask(magnet::function::Task::makeTask(addWindowFunc, window));

    //Spinlock waiting for the window to initialize
    while (!window->isReady()) { smallSleep(); }
  }

  inline void 
  shutdownCoil() 
  { 
    magnet::thread::ScopedLock lock(_coilLock);
    _runFlag = false; 
    _coilReadyFlag = false;
  }

  void waitForShutdown();

  inline bool isRunning() { return _runFlag; }

  //This mutex exists to stop coil killing itself while the main
  //program is accessing it. It is locked while windows are added or
  //destroyed and before shutting down.  It should be locked by the
  //main program while doing anything that requires a window to stay
  //alive. You should lock, then check isRunning() and abort the
  //action if it's not.
  magnet::thread::Mutex _coilLock;

private:
  CoilMaster();
  ~CoilMaster();
  
  void coilThreadEntryPoint();

  void smallSleep();

  volatile bool _runFlag; 
  volatile bool _coilReadyFlag;

  magnet::thread::Thread _coilThread;
  magnet::thread::TaskQueue _coilQueue;
   
  ///////////////////////////Glut GL render layer//////////////////////////

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
    
  ///////////////////////////GTK window layer/////////////////////////////
  static inline void addWindowFunc(CoilWindow* window)
  {
    magnet::thread::ScopedLock lock(CoilMaster::getInstance()._coilLock);

    window->init();
  }

  Gtk::Main _GTKit;

  bool GTKIldeFunc();
};

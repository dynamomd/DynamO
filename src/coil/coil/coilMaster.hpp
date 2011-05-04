/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include <magnet/thread/thread.hpp>
#include <magnet/thread/mutex.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/thread/refPtr.hpp>

#include <map>
#include <memory>

class CoilMaster {
public:
  //Only for window objects to call
  void  CallGlutCreateWindow(const char*, CoilWindow*);
  void  unregisterWindow(CoilWindow*);

  inline bool isRunning() { return _runFlag; }

  void addWindow(magnet::thread::RefPtr<CoilWindow>& window)
  {
    if (!isRunning()) M_throw() << "Coil is not running, cannot add a window";

    _coilQueue.queueTask(magnet::function::Task::makeTask(&CoilMaster::addWindowFunc, this, window));

    //Spinlock waiting for the window to initialize
    while (!window->isReady()) { smallSleep(); }
  }

  magnet::thread::TaskQueue& getTaskQueue() { return _coilQueue; }

  //This mutex exists to stop coil killing itself while the main
  //program is accessing it. It is locked while windows are added or
  //destroyed and before shutting down.  It should be locked by the
  //main program while doing anything that requires a window to stay
  //alive. You should lock, then check isRunning() and abort the
  //action if it's not.
  magnet::thread::Mutex _coilLock;

private:
  friend class CoilRegister;

  inline void 
  shutdownCoil() 
  { 
    magnet::thread::ScopedLock lock(_coilLock);
    _runFlag = false; 
    _coilReadyFlag = false;
  }

  void waitForShutdown();
  void bootRenderThread();

  CoilMaster();
  ~CoilMaster();
  
  void coilThreadEntryPoint();

  void smallSleep();

  volatile bool _runFlag; 
  volatile bool _coilReadyFlag;

  magnet::thread::Thread _coilThread;
  magnet::thread::TaskQueue _coilQueue;

  ///////////////////////////Glut GL render layer//////////////////////////

  std::map<int, magnet::thread::RefPtr<CoilWindow> > _viewPorts;

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
  inline void addWindowFunc(magnet::thread::RefPtr<CoilWindow>& window)
  {
    magnet::thread::ScopedLock lock(_coilLock);
    window->init();
    _viewPorts[window->GetWindowID()] = window;
  }



  std::auto_ptr<Gtk::Main> _GTKit;

  bool glutIdleTimeout();

  //This is a GTK timeout function to make sure any tasks for the coil
  //thread are performed. This is performed in a timer as it's too
  //expensive to do all the time (due to the lock)
  bool taskTimeout();
};


//!This class is like a connection register for the CoilMaster singleton.
class CoilRegister
{
public:
  inline CoilRegister() { increment(); }

  inline CoilRegister(const CoilRegister&) { increment(); }

  inline CoilRegister& operator=(const CoilRegister&) { return *this; }
  
  inline ~CoilRegister() { decrement(); }

  inline CoilMaster& getInstance() { return *_instance; }

private:
  friend class CoilMaster;
  friend class CLGLWindow;
  friend class RSphericalParticles;
  //! This instance is only for the CoilMaster and window classes to
  //! use, Everything else should use an instance of the register class
  //! to access coil!
  inline static 
  CoilMaster& getCoilInstance() { return *_instance; }
  
  inline void increment()
  {
    _mutex.lock();
    if (++_counter == 1)
      _instance = new CoilMaster;
    
    _mutex.unlock();
  }

  inline void decrement()
  {
    _mutex.lock();

    if (--_counter == 0)
      { delete _instance; _instance = NULL; }
    
    _mutex.unlock();
  }

  static CoilMaster* _instance;
  static magnet::thread::Mutex _mutex;
  static size_t _counter;
};

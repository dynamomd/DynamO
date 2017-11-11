/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <coil/coilWindow.hpp>
#include <memory>
#include <map>
#include <mutex>
#include <thread>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/exception.hpp>

namespace coil {
  class CoilMaster {
  public:
    //Only for window objects to call
    void  CallGlutCreateWindow(const char*, CoilWindow*);
    void  unregisterWindow(CoilWindow*);

    inline bool isRunning() { return _runFlag; }

    template<class T>
    void addWindow(std::shared_ptr<T>& window)
    {
      std::shared_ptr<CoilWindow> win 
	= std::static_pointer_cast<CoilWindow>(window);
      if (!isRunning()) M_throw() << "Coil is not running, cannot add a window";

      _coilQueue.queueTask(std::bind(&CoilMaster::addWindowFunc, this, win));

      if (parallel)
	//Spinlock waiting for the window to initialize
	while (parallel && (!window->isReady()))
	  {
	    std::this_thread::yield();
	    if (!isRunning())
	      M_throw() << "Coil failed to add the window as the main render thread died";
	  }
      else
	//Process the window addition now
	_coilQueue.drainQueue();
    }

    magnet::thread::TaskQueue& getTaskQueue() { return _coilQueue; }

    //This mutex exists to stop coil killing itself while the main
    //program is accessing it. It is locked while windows are added or
    //destroyed and before shutting down.  It should be locked by the
    //main program while doing anything that requires a window to stay
    //alive. You should lock, then check isRunning() and abort the
    //action if it's not.
    std::mutex _coilLock;

    static bool parallel;

    void init_tasks();
    bool main_loop_iter();
    void render_thread_entry_point();
    
  private:
    friend class CoilRegister;

    inline void 
    shutdownCoil() 
    { 
      std::lock_guard<std::mutex> lock(_coilLock);
      _runFlag = false; 
      _coilReadyFlag = false;
    }

    void waitForShutdown();

    CoilMaster();
    CoilMaster(bool);
    ~CoilMaster();
  
    void renderThreadShutdownTasks();

    volatile bool _runFlag; 
    volatile bool _coilReadyFlag;

    std::thread _coilThread;
    magnet::thread::TaskQueue _coilQueue;

    ///////////////////////////Glut GL render layer//////////////////////////

    std::map<int, std::shared_ptr<CoilWindow> > _viewPorts;

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
    inline void addWindowFunc(std::shared_ptr<CoilWindow> window)
    {
      std::lock_guard<std::mutex> lock(_coilLock);
      window->init();
      _viewPorts[window->GetWindowID()] = window;
    }

    std::unique_ptr<Gtk::Main> _GTKit;

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
    static std::mutex _mutex;
    static size_t _counter;
  };
}

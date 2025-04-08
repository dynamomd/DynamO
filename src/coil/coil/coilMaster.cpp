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
#include <GL/glew.h>
#include <gtkmm.h> // Must be included first!

#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#include <coil/coilMaster.hpp>
#include <coil/coilWindow.hpp>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <magnet/arg_share.hpp>

namespace coil {
CoilMaster *CoilRegister::_instance = NULL;
std::mutex CoilRegister::_mutex;
size_t CoilRegister::_counter = 0;
bool CoilMaster::parallel = true;

CoilMaster::CoilMaster() : _runFlag(true), _coilReadyFlag(false) {
  if (parallel) {
    _coilThread =
        std::thread(std::bind(&CoilMaster::render_thread_entry_point, this));

    // Spinlock waiting for the boot thread to come up
    while (!_coilReadyFlag) {
      std::this_thread::yield();
    }
  } else
    init_tasks();
}

CoilMaster::~CoilMaster() {
  // The thread must be shut down first before we can destroy the data members
  // of this class
  shutdownCoil();
  waitForShutdown();
}

void CoilMaster::CallBackDisplayFunc() {
  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackDisplayFunc();
}

void CoilMaster::CallBackCloseWindow() {
  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->deinit();
  CoilRegister::getCoilInstance()._viewPorts.erase(windowID);
}

void CoilMaster::CallBackKeyboardFunc(unsigned char key, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackKeyboardFunc(
      key, x, y);
}

void CoilMaster::CallBackKeyboardUpFunc(unsigned char key, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackKeyboardUpFunc(
      key, x, y);
}

void CoilMaster::CallBackMotionFunc(int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackMotionFunc(x,
                                                                           y);
}

void CoilMaster::CallBackMouseFunc(int button, int state, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackMouseFunc(
      button, state, x, y);
}

void CoilMaster::CallBackMouseWheelFunc(int button, int dir, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackMouseWheelFunc(
      button, dir, x, y);
}

void CoilMaster::CallBackPassiveMotionFunc(int x, int y) {
  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()
      ._viewPorts[windowID]
      ->CallBackPassiveMotionFunc(x, y);
}

void CoilMaster::CallBackReshapeFunc(int w, int h) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackReshapeFunc(w,
                                                                            h);
}

void CoilMaster::CallBackSpecialFunc(int key, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackSpecialFunc(
      key, x, y);
}

void CoilMaster::CallBackSpecialUpFunc(int key, int x, int y) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackSpecialUpFunc(
      key, x, y);
}

void CoilMaster::CallBackVisibilityFunc(int visible) {

  int windowID = glutGetWindow();

#ifdef DYNAMO_DEBUG
  if (!CoilRegister::getCoilInstance()._viewPorts.count(windowID))
    M_throw() << "Missing viewport!";
#endif

  CoilRegister::getCoilInstance()._viewPorts[windowID]->CallBackVisibilityFunc(
      visible);
}

void CoilMaster::CallGlutCreateWindow(const char *setTitle,
                                      CoilWindow *coilWindow) {

  // Open new window, record its windowID ,

  int windowID = glutCreateWindow(setTitle);

  coilWindow->SetWindowID(windowID);

  // Store the address of new window in global array
  // so CoilMaster can send events to proper callback functions.

  // Hand address of universal static callback functions to Glut.
  // This must be for each new window, even though the address are constant.
  glutDisplayFunc(CallBackDisplayFunc);

  // Idling is handled in coilMasters main loop
  // Timeout for render
  glutKeyboardFunc(CallBackKeyboardFunc);
  glutKeyboardUpFunc(CallBackKeyboardUpFunc);
  glutSpecialFunc(CallBackSpecialFunc);
  glutSpecialUpFunc(CallBackSpecialUpFunc);
  glutMouseFunc(CallBackMouseFunc);
  glutMouseWheelFunc(CallBackMouseWheelFunc);
  glutMotionFunc(CallBackMotionFunc);
  glutPassiveMotionFunc(CallBackPassiveMotionFunc);
  glutReshapeFunc(CallBackReshapeFunc);
  glutVisibilityFunc(CallBackVisibilityFunc);
  glutCloseFunc(CallBackCloseWindow);
}

void CoilMaster::unregisterWindow(CoilWindow *coilWindow) {
  int windowID = coilWindow->GetWindowID();
  _viewPorts.erase(windowID);
}

void CoilMaster::waitForShutdown() {
  if (parallel && _coilThread.joinable())
    _coilThread.join();
}

void CoilMaster::init_tasks() {
  _GTKit.reset(new Gtk::Main(magnet::ArgShare::getInstance().getArgc(),
                             magnet::ArgShare::getInstance().getArgv()));

  glutInit(&magnet::ArgShare::getInstance().getArgc(),
           magnet::ArgShare::getInstance().getArgv());

  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

  // Register the idle function
  // Glib::signal_idle().connect(sigc::mem_fun(this,
  // &CoilMaster::glutIdleTimeout));
  // Glib::signal_timeout().connect(sigc::mem_fun(this,
  // &CoilMaster::glutIdleTimeout), 30, Glib::PRIORITY_DEFAULT_IDLE);

  // This timeout is used to perform occasional tasks (that might block)
  Glib::signal_timeout().connect(sigc::mem_fun(this, &CoilMaster::taskTimeout),
                                 50, Glib::PRIORITY_DEFAULT_IDLE);

  _coilReadyFlag = true;
}

bool CoilMaster::main_loop_iter() {
  for (auto &win : _viewPorts)
    win.second->CallBackIdleFunc();

  glutMainLoopEvent();

  if (!parallel && _viewPorts.empty())
    return false;

  return _GTKit->iteration(false);
}

void CoilMaster::render_thread_entry_point() {
  try {
    init_tasks();
    while (main_loop_iter()) {
    };
  } catch (std::exception &except) {
    std::cerr << "\nRender thread caught an exception\n"
              << except.what() << "\n";

    shutdownCoil();
    renderThreadShutdownTasks();
  } catch (...) {
    std::cerr << "\nRender thread caught an unknown exception!\n";
    shutdownCoil();
    renderThreadShutdownTasks();
  }
}

void CoilMaster::renderThreadShutdownTasks() {
  // The task queue should not get any more tasks as
  // CoilRegister::getCoilInstance().isRunning() is false

  //! \todo{There is a race condition here if a window is added as coil is
  //! shutting down}
  {
    std::lock_guard<std::mutex> lock(_coilLock);

    // Now we must get glut to destroy all the windows

    while (!_viewPorts.empty()) {
      glutDestroyWindow(_viewPorts.begin()->first);
      // Run glutMainLoopEvent to let destroyed windows close
      glutMainLoopEvent();
      glutMainLoopEvent();
      glutMainLoopEvent();
      glutMainLoopEvent();
      glutMainLoopEvent();
      glutMainLoopEvent();
    }
  }

  Gtk::Main::quit();
}

bool CoilMaster::taskTimeout() {
  try {
    _coilQueue.drainQueue();

    if (!isRunning())
      renderThreadShutdownTasks();
  } catch (std::exception &except) {
    std::cerr << "\nCoil caught an exception while performing its tasks\n"
              << except.what() << "\n";
    shutdownCoil();
    renderThreadShutdownTasks();
  }

  return true;
}

// Somewhere for the vtable to live
CoilWindow::CoilWindow() : _readyFlag(false) {}
} // namespace coil

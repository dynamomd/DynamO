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

#include <magnet/arg_share.hpp>
#include <coil/coilMaster.hpp>
#include <coil/coilWindow.hpp>
#include <GL/freeglut_ext.h>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
                 
        
CoilMaster::CoilMaster():
  _runFlag(false),
  _renderReadyFlag(false),
  _windowReadyFlag(false),
  _GTKit(magnet::ArgShare::getInstance().getArgc(), magnet::ArgShare::getInstance().getArgv())
{
  glutInit(&magnet::ArgShare::getInstance().getArgc(), magnet::ArgShare::getInstance().getArgv());
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
}

CoilMaster::~CoilMaster(){
  shutdownCoil();
  waitForShutdown();
}
 
void CoilMaster::CallBackDisplayFunc(void){
   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackDisplayFunc();
}

void CoilMaster::CallBackCloseWindow()
{
  //Shutdown all windows, we can't recover from a window close in glut
  CoilMaster::getInstance().shutdownCoil();
}

void CoilMaster::CallBackKeyboardFunc(unsigned char key, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackKeyboardFunc(key, x, y);
}

void CoilMaster::CallBackKeyboardUpFunc(unsigned char key, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackKeyboardUpFunc(key, x, y);
}

void CoilMaster::CallBackMotionFunc(int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackMotionFunc(x, y);
}

void CoilMaster::CallBackMouseFunc(int button, int state, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackMouseFunc(button, state, x, y);
}

void CoilMaster::CallBackMouseWheelFunc(int button, int dir, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackMouseWheelFunc(button, dir, x, y);
}

void CoilMaster::CallBackPassiveMotionFunc(int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackPassiveMotionFunc(x, y);
}

void CoilMaster::CallBackReshapeFunc(int w, int h){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackReshapeFunc(w, h);
}

void CoilMaster::CallBackSpecialFunc(int key, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackSpecialFunc(key, x, y);
}   

void CoilMaster::CallBackSpecialUpFunc(int key, int x, int y){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackSpecialUpFunc(key, x, y);
}   

void CoilMaster::CallBackVisibilityFunc(int visible){

   int windowID = glutGetWindow();
   CoilMaster::getInstance()._viewPorts[windowID]->CallBackVisibilityFunc(visible);
}

void CoilMaster::CallGlutCreateWindow(const char * setTitle, CoilWindow * coilWindow){

   // Open new window, record its windowID , 

   int windowID = glutCreateWindow(setTitle);

   coilWindow->SetWindowID(windowID);

   // Store the address of new window in global array 
   // so CoilMaster can send events to propoer callback functions.

   CoilMaster::getInstance()._viewPorts[windowID] = coilWindow;
   CoilMaster::getInstance()._windows[coilWindow] = windowID;

   // Hand address of universal static callback functions to Glut.
   // This must be for each new window, even though the address are constant.
   glutDisplayFunc(CallBackDisplayFunc);
   
   //Idling is handled in coilMasters main loop
   glutIdleFunc(NULL);
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

void CoilMaster::CallGlutDestroyWindow(CoilWindow * coilWindow)
{
  glutDestroyWindow(_windows[coilWindow]);
}

void CoilMaster::renderThreadEntryPoint()
{
  try {
    _renderReadyFlag = true;

    while (CoilMaster::getInstance()._runFlag)
	{
	  //Don't just burn cycles if we have no windows!
	  if (CoilMaster::getInstance()._viewPorts.empty())
	    smallSleep();
	  else
	    glutMainLoopEvent();
	  
	  
	  for (std::map<int,CoilWindow*>::iterator iPtr = CoilMaster::getInstance()._viewPorts.begin();
	       iPtr != CoilMaster::getInstance()._viewPorts.end(); ++iPtr)
	    {
	      glutSetWindow(iPtr->first);
	      iPtr->second->CallBackIdleFunc();
	    }

	  //Now we drain our task queue
	  _renderQueue.drainQueue();
	}


    {
      magnet::thread::ScopedLock lock(_coilLock);
      //If we reach here, we must just delete all windows that we own to
      //free up the memory. There is a danger that the renderQueue has
      //some non-zero size, with a window left to initialize.
      //! \todo{Fix the race condition on adding windows and deleting them}
      for (std::map<int,CoilWindow*>::iterator iPtr = CoilMaster::getInstance()._viewPorts.begin();
	   iPtr != CoilMaster::getInstance()._viewPorts.end(); ++iPtr)
	delete iPtr->second;
      
      _viewPorts.clear();
      _windows.clear();
    }

    //Run glutMainLoopEvent to let destroyed windows close
    glutMainLoopEvent();
  } catch (std::exception& except)
    {
      std::cerr << "\nRender thread caught an exception\n"
		<< except.what()
		<< "\n As we're in a thread we can only exit(1)!";
      std::exit(1);
    } catch (...)
    {
      std::cerr << "\nRender thread caught an unknown exception!\n"
		<< "\n As we're in a thread we can only exit(1)!";
      std::exit(1);
    }
}

void CoilMaster::smallSleep()
{
  timespec sleeptime;
  sleeptime.tv_sec = 0;
  sleeptime.tv_nsec = 100000000;
  nanosleep(&sleeptime, NULL);
}

void 
CoilMaster::bootCoil()
{
  if (_runFlag)
    throw std::runtime_error("Renderer is already running!");

  _runFlag = true;
  
  _renderThread = magnet::thread::Thread(magnet::function::Task::makeTask(&CoilMaster::renderThreadEntryPoint, this));
  _windowThread = magnet::thread::Thread(magnet::function::Task::makeTask(&CoilMaster::windowThreadEntryPoint, this));

  //Spinlock waiting for the boot thread to come up
  while (!_renderReadyFlag || !_windowReadyFlag) { smallSleep(); }
}

void 
CoilMaster::waitForShutdown()
{
  if (_renderThread.validTask()) _renderThread.join();
  if (_windowThread.validTask()) _windowThread.join();
}

//The glade xml file is "linked" into a binary file and stuffed in the executable, these are the symbols to its data
extern const char _binary_src_coil_coil_gui_gladexml_start[];
extern const char _binary_src_coil_coil_gui_gladexml_end[];

void CoilMaster::windowThreadEntryPoint()
{
  try {
    //Build the windows from the glade data
    {
      Glib::ustring glade_data(reinterpret_cast<const char *>(_binary_src_coil_coil_gui_gladexml_start), 
			       _binary_src_coil_coil_gui_gladexml_end
			       -_binary_src_coil_coil_gui_gladexml_start);
      
      _refXml = Gtk::Builder::create_from_string(glade_data);
    }
    
    
    //Register the idle function
    Glib::signal_timeout().connect(sigc::mem_fun(this, &CoilMaster::GTKIldeFunc), 100);
    
    {//Now setup some callback functions
      Gtk::Window* controlwindow;
      _refXml->get_widget("controlWindow", controlwindow);
      
      //Setup the on_close button
      controlwindow->signal_delete_event().connect(sigc::mem_fun(this, &CoilMaster::onControlWindowDelete));
    }

    _windowReadyFlag = true;
    _GTKit.run();
  } catch (std::exception& except)
    {
      std::cerr << "\nRender thread caught an exception\n"
		<< except.what()
		<< "\n As we're in a thread we can only exit(1)!";
      std::exit(1);
    } catch (...)
    {
      std::cerr << "\nRender thread caught an unknown exception!\n"
		<< "\n As we're in a thread we can only exit(1)!";
      std::exit(1);
    }
}

bool CoilMaster::onControlWindowDelete(GdkEventAny * pEvent)
{
  shutdownCoil();

  return false;
}

bool CoilMaster::GTKIldeFunc()
{
  if (!CoilMaster::getInstance().isRunning()) 
    {
      Gtk::Main::quit();
      Gtk::Window* controlwindow;
      _refXml->get_widget("controlWindow", controlwindow);
      
      //Setup the on_close button
      controlwindow->hide();
    }

  return true;
}

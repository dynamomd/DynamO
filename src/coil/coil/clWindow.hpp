/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include <vector>
#include <gtkmm.h>
#include <coil/coilMaster.hpp>
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <coil/Maths/Maths.h>
#include <magnet/math/vector.hpp>
#include <magnet/static_assert.hpp>

#include <magnet/GL/light.hpp>
#include <magnet/CL/GLBuffer.hpp>
#include <magnet/GL/shadowShader.hpp>
#include <magnet/GL/shadowFBO.hpp>
#include <magnet/GL/viewPort.hpp>
#include <magnet/GL/multisampledFBO.hpp>
#include <magnet/GL/nrmlShader.hpp>
#include <magnet/thread/refPtr.hpp>

#include <coil/filters/filter.hpp>
#include <coil/RenderObj/RenderObj.hpp>
#include <memory>

#ifdef COIL_wiimote
# include <coil/extcode/wiiheadtracking.hpp>
#endif 

class CLGLWindow : public CoilWindow
{
public:
  CLGLWindow(int setWidth, int setHeight,
	     int setInitPositionX, int setInitPositionY,
	     std::string title, double updateIntervalValue,
	     bool dynamo = false);

  ~CLGLWindow();
  
  virtual void CallBackDisplayFunc();
  virtual bool CallBackIdleFunc();

  void CallBackReshapeFunc(int w, int h);    

  const std::string& getWindowTitle() const { return windowTitle; }
  void setWindowtitle(const std::string& newtitle);
  
  void addRenderObj(const magnet::thread::RefPtr<RenderObj>& nObj)
  { RenderObjects.push_back(nObj); }

  inline volatile const int& getLastFrameTime() const { return _lastFrameTime; }

  magnet::CL::CLGLState& getCLState() { return *_CLState; }

  void init();
  void deinit();

  magnet::thread::Mutex& getDestroyLock() { return _destroyLock; }

  bool simupdateTick();
  
  const double& getUpdateInterval() {return _updateIntervalValue; }

  void flagNewData() { _newData = true; }

  void setUpdateRateUnitToSteps(size_t defaultsteps = 100);

  void setSimStatus1(std::string);
  void setSimStatus2(std::string);

  magnet::thread::RefPtr<magnet::thread::TaskQueue>& getQueue() { return  _systemQueue; }
protected:
  void setLabelText(Gtk::Label*, std::string);

  magnet::GL::shadowShader _shadowShader;
  magnet::GL::shadowFBO _shadowFBO;

  //Primary render target
  std::auto_ptr<magnet::GL::FBO> _renderTarget;

  //Frame buffers to flip flop filters between
  magnet::GL::FBO _filterTarget1;
  magnet::GL::FBO _filterTarget2;

  magnet::GL::NormalShader _nrmlShader;
  magnet::GL::FBO _normalsFBO;

  size_t _height, _width;
  int _windowX, _windowY;

#ifdef COIL_wiimote
  TrackWiimote _wiiMoteTracker;
#endif

  std::vector<magnet::thread::RefPtr<RenderObj> > RenderObjects;

  void CallBackSpecialUpFunc(int key, int x, int y) {}
  void CallBackSpecialFunc(int key, int x, int y) {} 
  void CallBackKeyboardFunc(unsigned char key, int x, int y);
  void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
  void CallBackMouseWheelFunc(int button, int dir, int x, int y);
  void CallBackMouseFunc(int button, int state, int x, int y);
  void CallBackMotionFunc(int x, int y);

  void performPicking(int x, int y);



private:
  //Task queue for the simulation thread
  magnet::thread::RefPtr<magnet::thread::TaskQueue> _systemQueue;
  double _updateIntervalValue;

  magnet::thread::Mutex _destroyLock;

  void CameraSetup();

  magnet::thread::RefPtr<magnet::CL::CLGLState> _CLState;

  
  virtual void initOpenGL();
  virtual void initOpenCL();

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
  size_t _frameCounter, _updateCounter;

  volatile int _lastFrameTime;
  int _FPStime; 
  int _lastUpdateTime;
  int _frameRenderTime;

  sigc::connection _renderTimeout;

  std::auto_ptr<magnet::GL::viewPort> _viewPortInfo;
    
  bool keyStates[256];

  float _mouseSensitivity; 
  float _moveSensitivity;
 
  int _oldMouseX, _oldMouseY;
  int _specialKeys;

  bool _shaderPipeline;
  bool _shadowMapping;
  GLfloat _shadowIntensity;
  volatile bool _simrun;
  volatile bool _simframelock;
  bool _snapshot;
  bool _record;
  bool _showLight; 
  bool _PNGFileFormat;
  bool _fpsLimit;
  int  _fpsLimitValue;
  bool _filterEnable;
  bool _pickingEnabled;

  size_t _snapshot_counter;

  std::auto_ptr<magnet::GL::lightInfo> _light0;

  /////////GTK members
  virtual void initGTK();

  bool GTKTick();
  Glib::RefPtr<Gtk::Builder> _refXml;
  Gtk::Window* controlwindow;
  sigc::connection _timeout_connection;

  struct FilterModelColumnsType : Gtk::TreeModelColumnRecord
  {
    FilterModelColumnsType()
    { add(m_active); add(m_name); add(m_filter_ptr);}
    
    Gtk::TreeModelColumn<bool> m_active;
    Gtk::TreeModelColumn<Glib::ustring> m_name;
    Gtk::TreeModelColumn<void*> m_filter_ptr;
  };

  struct RenderObjModelColumnsType : Gtk::TreeModelColumnRecord
  {
    RenderObjModelColumnsType()
    { add(m_name); add(m_visible); add(m_id);}
    
    Gtk::TreeModelColumn<Glib::ustring> m_name;
    Gtk::TreeModelColumn<bool> m_visible;
    Gtk::TreeModelColumn<size_t> m_id;
  };

  std::auto_ptr<FilterModelColumnsType> _filterModelColumns;
  std::auto_ptr<RenderObjModelColumnsType> _renderObjModelColumns;

  Glib::RefPtr<Gtk::ListStore> _filterStore;
  Glib::RefPtr<Gtk::ListStore> _renderObjStore;

  Gtk::TreeView* _filterView;
  Gtk::TreeView* _renderObjView;

  //Callback for enabling/disabling the shader pipeline
  void pipelineEnableCallback();

  //Filter control callbacks
  void filterUpCallback();
  void filterDownCallback();
  void filterDeleteCallback();
  void filterAddCallback();
  void filterSelectCallback();
  void filterClearCallback();
  void filterActiveCallback();
  
  //Render Object Contol Callbacks
  void rebuildRenderView();
  void visibleRObjCallback();
  void editRObjCallback();
  void deleteRObjCallback();
  void addRObjCallback();
  void selectRObjCallback();


  //Other callbacks
  void multisampleEnableCallback();
  void shadowEnableCallback();

  void simRunControlCallback();
  void simFramelockControlCallback();
  void snapshotCallback();
  void recordCallback();
  void axisShowCallback();
  void lightShowCallback();
  void lightPlaceCallback();
  void shadowIntensityCallback(double);
  void snapshotFileFormatCallback();
  void FPSLimitCallback();
  void aboutCallback();
  void renderNormalsCallback();
  void renderModeCallback();

  void runCallback(); 

  //Generic Value update callback
  void guiUpdateCallback();

  //Dynamo specifc stuff
public:
  inline bool dynamoParticleSync() const { return _particleSync; }

private:
  bool _dynamo;
  bool _particleSync;
  volatile bool _newData;

};

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
#include <magnet/math/vector.hpp>
#include <magnet/static_assert.hpp>
#include <magnet/GL/light.hpp>
#include <magnet/GL/shader/render.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/multisampledFBO.hpp>
#include <magnet/GL/shader/normal.hpp>
#include <magnet/GL/shader/vsm.hpp>
#include <magnet/GL/shader/simple_render.hpp>
#include <coil/filters/filter.hpp>
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/coilMaster.hpp>
#include <boost/signals2.hpp>
#include <vector>
#include <memory>

namespace coil {
  class CLGLWindow : public CoilWindow
  {
  public:
    CLGLWindow(std::string title, double updateIntervalValue,
	       bool dynamo = false);

    ~CLGLWindow();
  
    virtual void CallBackDisplayFunc();
    virtual bool CallBackIdleFunc();

    void CallBackReshapeFunc(int w, int h);    

    const std::string& getWindowTitle() const { return windowTitle; }
    void setWindowtitle(const std::string& newtitle);
  
    void addRenderObj(const std::tr1::shared_ptr<RenderObj>& nObj)
    { _renderObjsTree._renderObjects.push_back(nObj); }

    inline volatile const int& getLastFrameTime() const { return _lastFrameTime; }

    void init();
    void deinit();

    void simupdateTick(double t);
  
    const double& getUpdateInterval() {return _updateIntervalValue; }

    void setUpdateRateUnitToSteps(size_t defaultsteps = 100);

    void setSimStatus1(std::string);
    void setSimStatus2(std::string);

    std::tr1::shared_ptr<magnet::thread::TaskQueue>& getQueue() { return  _systemQueue; }

    magnet::GL::Context& getGLContext()
    {
      if (_glContext == NULL)
	M_throw() << "GL context not yet set!";
      return *_glContext;
    }

    boost::signals2::signal<void ()>& signal_data_update() { return _updateDataSignal; }

  protected:
    boost::signals2::signal<void ()> _updateDataSignal;
    
    void setLabelText(Gtk::Label*, std::string);

    magnet::GL::shader::RenderShader _renderShader;
    magnet::GL::shader::VSMShader _VSMShader;
    magnet::GL::shader::SimpleRenderShader _simpleRenderShader;

    //Primary render target
    std::auto_ptr<magnet::GL::FBO> _renderTarget;

    //Frame buffers to flip flop filters between
    magnet::GL::FBO _filterTarget1;
    magnet::GL::FBO _filterTarget2;

    magnet::GL::shader::NormalShader _nrmlShader;
    magnet::GL::FBO _normalsFBO;

    void CallBackSpecialUpFunc(int key, int x, int y) {}
    void CallBackSpecialFunc(int key, int x, int y) {} 
    void CallBackKeyboardFunc(unsigned char key, int x, int y);
    void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
    void CallBackMouseWheelFunc(int button, int dir, int x, int y);
    void CallBackMouseFunc(int button, int state, int x, int y);
    void CallBackMotionFunc(int x, int y);

    void performPicking(int x, int y);

    //Task queue for the simulation thread
    std::tr1::shared_ptr<magnet::thread::TaskQueue> _systemQueue;
    double _updateIntervalValue;
    size_t _consoleID;

    magnet::GL::Context* _glContext;

    magnet::thread::Mutex _destroyLock;

    void CameraSetup();

    void drawScene(magnet::GL::FBO&, magnet::GL::Camera&);

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

    magnet::GL::Camera _camera;
    
    bool keyStates[256];

    float _mouseSensitivity; 
    float _moveSensitivity;
 
    int _oldMouseX, _oldMouseY;
    int _specialKeys;

    bool _shadowMapping;
    GLfloat _shadowIntensity;
    volatile bool _simrun;
    volatile bool _simframelock;
    bool _snapshot;
    bool _record;
    bool _PNGFileFormat;
    bool _fpsLimit;
    int  _fpsLimitValue;
    bool _filterEnable;
    bool _analygraphMode;

    size_t _snapshot_counter;

    magnet::GL::Light _light0;

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

    std::auto_ptr<FilterModelColumnsType> _filterModelColumns;

    Glib::RefPtr<Gtk::ListStore> _filterStore;
    Gtk::TreeView* _filterView;

    RenderObjectsGtkTreeView _renderObjsTree;

    //Filter control callbacks
    void filterUpCallback();
    void filterDownCallback();
    void filterDeleteCallback();
    void filterAddCallback();
    void filterSelectCallback();
    void filterClearCallback();
    void filterActiveCallback();
  
    //Render Object Contol Callbacks
    void selectRObjCallback();

    //Wii Remote callbacks
    void wiiMoteConnect(); 
    bool wiiMoteIRExposeEvent(GdkEventExpose*);

    void HeadReset();

    //Other callbacks
    void multisampleEnableCallback();
    void shadowEnableCallback();

    void simRunControlCallback();
    void simFramelockControlCallback();
    void snapshotCallback();
    void recordCallback();
    void axisShowCallback();
    void lightPlaceCallback();
    void shadowIntensityCallback(double);
    void snapshotFileFormatCallback();
    void FPSLimitCallback();
    void aboutCallback();
    void renderNormalsCallback();

    void runCallback(); 

    //Generic Value update callback
    void guiUpdateCallback();

    //Dynamo specifc stuff
  public:
    inline bool dynamoParticleSync() const { return _particleSync; }

  protected:
    bool _dynamo;
    bool _particleSync;
    volatile bool _newData;

  };
}

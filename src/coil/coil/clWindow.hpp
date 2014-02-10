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
#include <magnet/math/vector.hpp>
#include <magnet/GL/shader/lightShader.hpp>
#include <magnet/GL/shader/ambientLight.hpp>
#include <magnet/GL/shader/luminance.hpp>
#include <magnet/GL/shader/blur.hpp>
#include <magnet/GL/shader/toneMap.hpp>
#include <magnet/GL/shader/depthResolver.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/multisampledFBO.hpp>
#include <magnet/GL/shader/copy.hpp>
#include <magnet/GL/shader/downsampler.hpp>
#include <magnet/GL/objects/cairo.hpp>
#include <magnet/function/delegate.hpp>
#include <coil/filters/filter.hpp>
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/coilMaster.hpp>
#include <vector>
#include <memory>

namespace magnet {
  namespace image {
    class VideoEncoderFFMPEG;
  }
}

namespace coil {
  class CLGLWindow : public CoilWindow
  {
  public:
    CLGLWindow(std::string title, double updateIntervalValue, bool dynamo = false);

    ~CLGLWindow();
  
    virtual void CallBackDisplayFunc();
    virtual bool CallBackIdleFunc();

    void resizeRender(int w, int h);
    void CallBackReshapeFunc(int w, int h);    

    const std::string& getWindowTitle() const { return windowTitle; }
    void setWindowtitle(const std::string& newtitle);
  
    void addRenderObj(const std::shared_ptr<RenderObj>& nObj)
    { _renderObjsTree._renderObjects.push_back(nObj); }

    inline volatile const int& getLastFrameTime() const { return _lastFrameTime; }

    void init();
    void deinit();

    void simupdateTick(double t);
  
    const double& getUpdateInterval() {return _updateIntervalValue; }

    void setUpdateRateUnitToSteps(size_t defaultsteps = 100);

    void setSimStatus1(std::string);
    void setSimStatus2(std::string);

    std::shared_ptr<magnet::thread::TaskQueue>& getQueue() { return  _systemQueue; }

    magnet::GL::Context::ContextPtr& getGLContext()
    {
      if (!_glContext)
	M_throw() << "GL context not yet set!";
      return _glContext;
    }

    magnet::Signal<void()> _updateDataSignal;

    void autoscaleView();

    magnet::GL::Camera& getCamera() { return _camera; }
  protected:
    CLGLWindow(const CLGLWindow&);
    
    void setLabelText(Gtk::Label*, std::string);

    magnet::GL::shader::PointLightShader _pointLightShader;
    magnet::GL::shader::ShadowLightShader _shadowLightShader;
    magnet::GL::shader::AmbientLightShader _ambientLightShader;
    magnet::GL::shader::LuminanceShader _luminanceShader;
    magnet::GL::shader::LuminanceMipMapShader _luminanceMipMapShader;
    magnet::GL::shader::ToneMapShader _toneMapShader;
    magnet::GL::shader::DownsamplerShader _downsampleShader;
    magnet::GL::shader::SeperableGaussian _blurShader;
    magnet::GL::shader::DepthResolverShader _depthResolverShader;
    magnet::GL::shader::CopyShader _copyShader;

    //Primary render target, or the render target for the left eye.
    magnet::GL::FBO _renderTarget;
    magnet::GL::FBO _Gbuffer;
    magnet::GL::FBO _shadowbuffer;
    magnet::GL::FBO _hdrBuffer;
    magnet::GL::FBO _luminanceBuffer1;
    magnet::GL::FBO _luminanceBuffer2;

    //Blur Targets
    magnet::GL::FBO _blurTarget1;
    magnet::GL::FBO _blurTarget2;


    //Frame buffers to flip flop between
    magnet::GL::FBO _filterTarget1;
    magnet::GL::FBO _filterTarget2;

    //For object selection
    /*! \brief If valid, the render object which is currently
        selected. */
    std::shared_ptr<RenderObj> _selectedObject;
    /*! \brief If \ref _selectedObject is valid, this holds the id of
      the object within \ref _selectedObject which has been
      selected. */
    uint32_t _selectedObjectID;
    magnet::GL::objects::CairoSurface _cairo_screen;

    void CallBackSpecialUpFunc(int key, int x, int y) {}
    void CallBackSpecialFunc(int key, int x, int y) {} 
    void CallBackKeyboardFunc(unsigned char key, int x, int y);
    void CallBackKeyboardUpFunc(unsigned char key, int x, int y);
    void CallBackMouseWheelFunc(int button, int dir, int x, int y);
    void CallBackMouseFunc(int button, int state, int x, int y);
    void CallBackMotionFunc(int x, int y);

    void performPicking(int x, int y);

    //Task queue for the simulation thread
    std::shared_ptr<magnet::thread::TaskQueue> _systemQueue;
    double _updateIntervalValue;
    size_t _consoleID;

    magnet::GL::Context::ContextPtr _glContext;

    std::mutex _destroyLock;

    void CameraSetup();

    void drawScene(magnet::GL::Camera&);

    enum KeyStateType
      {
	DEFAULT = 0,
	LEFTMOUSE = 1,
	RIGHTMOUSE = 2,
	MIDDLEMOUSE = 4
      };

    int _mouseState;
  
    std::string windowTitle;
    bool FPSmode;
    size_t _frameCounter, _updateCounter;

    volatile int _lastFrameTime;
    int _FPStime; 
    int _lastUpdateTime;
    int _frameRenderTime;

    sigc::connection _renderTimeout;

    /*! \brief This camera must be a static member of the windows as
        other threads might queue tasks around it.
     */
    magnet::GL::Camera _camera;
    
    bool keyStates[256];

    float _mouseSensitivity; 
    float _moveSensitivity;
 
    int _oldMouseX, _oldMouseY;
    int _specialKeys;

    volatile bool _simrun;
    volatile bool _simframelock;
    bool _snapshot;
    bool _record;
    bool _PNGFileFormat;
    bool _fpsLimit;
    int  _fpsLimitValue;
    bool _filterEnable;
    bool _stereoMode;
    double _ambientIntensity;
    std::array<GLfloat, 3> _backColor;
    float _sceneKey;

    bool _bloomEnable;
    float _bloomCompression;
    float _bloomCutoff;
    float _bloomSaturation;

    size_t _snapshot_counter;
    size_t _video_counter;
    size_t _samples;

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

    std::unique_ptr<FilterModelColumnsType> _filterModelColumns;

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
    void simRunControlCallback();
    void simFramelockControlCallback();
    void snapshotCallback();
    void recordCallback();
    void axisShowCallback();
    void FPSLimitCallback();
    void aboutCallback();
    void renderNormalsCallback();

    void AAsamplechangeCallback();

    void addLightCallback();

    void addFunctionCallback();

    void cameraModeCallback();
    void runCallback(); 

    void rescaleCameraCallback();

    //Generic Value update callback
    void guiUpdateCallback();

    //Dynamo specifc stuff
  public:
    inline bool dynamoParticleSync() const { return _particleSync; }

  protected:
    bool _dynamo;
    bool _particleSync;
    volatile bool _newData;
    magnet::math::Vector _cameraFocus;
    enum CAM_MODE { ROTATE_WORLD, ROTATE_POINT, ROTATE_CAMERA};
    CAM_MODE _cameraMode;
    std::unique_ptr<Gtk::ComboBoxText> _aasamples;
#ifdef MAGNET_FFMPEG_SUPPORT
    std::unique_ptr<magnet::image::VideoEncoderFFMPEG> _encoder;
#endif
  };
}

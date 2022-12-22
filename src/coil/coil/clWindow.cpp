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

#include <coil/clWindow.hpp>
#include <coil/RenderObj/Surface.hpp>
#include <coil/RenderObj/console.hpp>
#include <coil/RenderObj/Volume.hpp>
#include <coil/RenderObj/Light.hpp>
#include <magnet/GL/context.hpp>
#include <magnet/image/videoEncoderFFMPEG.hpp>
#include <magnet/image/PNG.hpp>
#include <magnet/image/bitmap.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <gtkmm/volumebutton.h>
#include <stator/xml.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <locale>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef COIL_wiimote
# include <magnet/wiiheadtracking.hpp>
#endif 

#include <coil/images/images.hpp>

//The glade xml file is "linked" into a binary file and stuffed in the
//executable, these are the symbols to its data
extern const std::string clwingtk;

namespace {
  void resizeGlutWindow(int width, int height)
  {
    glutReshapeWindow(width, height);
    glutPostRedisplay();
  }
}

namespace coil {
  CLGLWindow::CLGLWindow(std::string title,
			 double updateIntervalValue,
			 bool dynamo
			 ):
    _selectedObjectID(0),
    _systemQueue(new magnet::thread::TaskQueue),
    _updateIntervalValue(updateIntervalValue),
    _mouseState(DEFAULT),
    windowTitle(title),
    _frameCounter(0),
    _updateCounter(0),
    _mouseSensitivity(0.3),
    _moveSensitivity(0.01),
    _specialKeys(0),
    _simrun(false),
    _simframelock(false),
    _snapshot(false),
    _record(false),
    _PNGFileFormat(true),
    _fpsLimit(false),
    _fpsLimitValue(25),
    _filterEnable(true),
    _stereoMode(false),
    _openVRMode(false),
    _ambientIntensity(0.001),
    _snapshot_counter(0),
    _video_counter(0),
    _samples(1),
    _dynamo(dynamo),
    _cameraMode(ROTATE_WORLD)
  {
    for (size_t i(0); i < 3; ++i)
      _backColor[i] = 0.0;

    for (size_t i(0); i < 256; ++i) keyStates[i] = false;
  }

  CLGLWindow::~CLGLWindow() {}

  bool
  CLGLWindow::CallBackIdleFunc()
  {
    try 
      {
	glutSetWindow(windowID);
	CallBackDisplayFunc();
      }
    catch (std::exception& except)
      {
	std::cerr << "\n Window render caught a std::exception\n"
		  << except.what();
	std::exit(1);
      }  
    catch (...)
      {
	std::cerr << "\nRender thread caught an unknown exception!\n";
	std::exit(1);
      }

    return true;
  }

  namespace {
    void setIcon(Glib::RefPtr<Gtk::Builder>& builder, std::string imgname, Glib::RefPtr<Gdk::Pixbuf> img)
    { 
      Gtk::Image* icon; 
      builder->get_widget(imgname, icon); 
      icon->set(img); 
    }
  }

  void
  CLGLWindow::initGTK()
  {
    _filterModelColumns.reset(new FilterModelColumnsType);

    try {
      _refXml = Gtk::Builder::create_from_string(clwingtk);
    } catch (std::exception& err)
      {
	M_throw() << "Failed to load the interface design into Gtk::Builder\n"
		  << err.what();
      }
	
    /////////Timeout for FPS and UPS calculation
    _timeout_connection = Glib::signal_timeout().connect_seconds(sigc::mem_fun(this, &CLGLWindow::GTKTick), 1);

    ////////Store the control window
    _refXml->get_widget("controlWindow", controlwindow);
  
    ////////Setup the window icon
    controlwindow->set_icon(coil::images::coilicon());
    setIcon(_refXml, "CamPlusXimg", coil::images::camplusx());
    setIcon(_refXml, "CamPlusYimg", coil::images::camplusy());
    setIcon(_refXml, "CamPlusZimg", coil::images::camplusz());
    setIcon(_refXml, "CamNegXimg",  coil::images::camnegx());
    setIcon(_refXml, "CamNegYimg",  coil::images::camnegy());
    setIcon(_refXml, "CamNegZimg",  coil::images::camnegz());
    setIcon(_refXml, "CamRescaleimg", coil::images::camrescale());
    setIcon(_refXml, "aboutSplashImage", coil::images::coilsplash());
    setIcon(_refXml, "addLightImage", coil::images::addLight_Icon());
    setIcon(_refXml, "addFunctionImage", coil::images::addFunction_Icon());

    {
      Gtk::Button* button;

      _refXml->get_widget("CamPlusXbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{1,0,0}));

      _refXml->get_widget("CamPlusYbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{0,1,0}));

      _refXml->get_widget("CamPlusZbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{0,0,1}));


      _refXml->get_widget("CamNegXbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{-1,0,0}));

      _refXml->get_widget("CamNegYbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{0,-1,0}));

      _refXml->get_widget("CamNegZbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::CameraHeadTracking::setViewAxis), magnet::math::Vector{0,0,-1}));
    }

    for (const auto& b : std::vector<std::pair<std::string, void(CLGLWindow::*)()> >{
	{"CamMode", &CLGLWindow::cameraModeCallback},
	{"LoadDataButton", &CLGLWindow::LoadDataCallback},
	{"SaveDataButton", &CLGLWindow::SaveDataCallback},
	{"addLightButton", &CLGLWindow::addLightCallback},
	{"addFunctionButton", &CLGLWindow::addFunctionCallback},
	{"SimSnapshot", &CLGLWindow::snapshotCallback},
	{"filterUp", &CLGLWindow::filterUpCallback},
	{"filterDown", &CLGLWindow::filterDownCallback},
	{"filterDelete", &CLGLWindow::filterDeleteCallback},
	{"filterAdd",&CLGLWindow::filterAddCallback},
	{"filterClear", &CLGLWindow::filterClearCallback},
	{"HeadTrackReset", &CLGLWindow::HeadReset}
      })
      {
	Gtk::Button* button;
	_refXml->get_widget(b.first, button);
	button->signal_clicked().connect(sigc::mem_fun(*this, b.second));
      }
    
    ///////Create a top menu
    {
      Glib::RefPtr<Gtk::UIManager> m_refUIManager;
      Glib::RefPtr<Gtk::ActionGroup> m_refActionGroup = Gtk::ActionGroup::create();

      m_refActionGroup->add(Gtk::Action::create("ViewMenu", "View"));
      m_refActionGroup->add(Gtk::Action::create("HelpMenu", "Help"));
      m_refActionGroup->add(Gtk::Action::create("AboutItem", "About"), sigc::mem_fun(this, &CLGLWindow::aboutCallback));
      m_refActionGroup->add(Gtk::Action::create("ResizeMenu", "Resize render window"));
      m_refActionGroup->add(Gtk::Action::create("800by600", " 800x600"), sigc::bind(&resizeGlutWindow, 800, 600));
      m_refActionGroup->add(Gtk::Action::create("1024by768", " 1024x768"), sigc::bind(&resizeGlutWindow, 1024, 768));
      m_refActionGroup->add(Gtk::Action::create("1280by720", " 1280x720"), sigc::bind(&resizeGlutWindow, 1280, 720));
      m_refActionGroup->add(Gtk::Action::create("1920by1080", " 1920x1080"), sigc::bind(&resizeGlutWindow, 1920, 1080));


      m_refUIManager = Gtk::UIManager::create();
      m_refUIManager->insert_action_group(m_refActionGroup);

      Glib::ustring ui_info = 
        "<ui>"
        "  <menubar name='MenuBar'>"
        "    <menu action='ViewMenu'>"
        "      <menu action='ResizeMenu'>"
	"        <menuitem action='800by600' />"
	"        <menuitem action='1024by768' />"
	"        <menuitem action='1280by720' />"
	"        <menuitem action='1920by1080' />"
        "      </menu>"
        "    </menu>"
	"    <menu action='HelpMenu'>"
        "      <menuitem action='AboutItem' />"
	"    </menu>"
        "  </menubar>"
        "</ui>";

      try
	{
	  m_refUIManager->add_ui_from_string(ui_info);
	}
      catch(const Glib::Error& ex)
	{
	  std::cerr << "building menus failed: " <<  ex.what();
	}

      Gtk::Widget* pMenubar = m_refUIManager->get_widget("/MenuBar");

      Gtk::VBox* pagebox;
      _refXml->get_widget("PageBox", pagebox);

      pagebox->pack_start(*pMenubar, Gtk::PACK_SHRINK);
      pagebox->reorder_child(*pMenubar, 0);
      pagebox->show_all_children();
    }  

    {////////Simulation run control
      Gtk::ToggleButton* togButton;
      _refXml->get_widget("SimRunButton", togButton);

      togButton->signal_toggled()
	.connect(sigc::mem_fun(this, &CLGLWindow::runCallback));
    }

    {///////Frame lock control
      Gtk::ToggleButton* framelockButton;
      _refXml->get_widget("SimLockButton", framelockButton);
      framelockButton->signal_toggled()
	.connect(sigc::mem_fun(this, &CLGLWindow::simFramelockControlCallback));
    }

    {
      namespace fs = boost::filesystem;
      Gtk::FileChooserButton* fileChooser;
      _refXml->get_widget("snapshotDirectory", fileChooser);
      fileChooser->set_filename(".");
    }

    {///////Recording button
      Gtk::ToggleButton* recordButton;
      _refXml->get_widget("SimRecordButton", recordButton);

#ifndef MAGNET_FFMPEG_SUPPORT
      recordButton->set_sensitive(false);
#endif

      recordButton->signal_toggled()
	.connect(sigc::mem_fun(this, &CLGLWindow::recordCallback));
    }
  
    {///////Control the update rate from the simulation
      Gtk::Entry* updateFreq;
      _refXml->get_widget("updateFreq", updateFreq);
      updateFreq->set_text
	(boost::lexical_cast<std::string>(_updateIntervalValue));
      updateFreq->signal_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
    }

    {///////FPS lock
      Gtk::ToggleButton* fpslockButton;
      _refXml->get_widget("FPSLimit", fpslockButton);
      fpslockButton->set_active(_fpsLimit);
      fpslockButton->signal_toggled()
	.connect(sigc::mem_fun(this, &CLGLWindow::FPSLimitCallback));
    }

    {///////FPS lock value
      Gtk::SpinButton* fpsButton;
      _refXml->get_widget("FPSLimitVal", fpsButton);
      fpsButton->set_value(_fpsLimitValue);
      fpsButton->signal_value_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::FPSLimitCallback));
    }


    ///////////////////////Render Pipeline//////////////////////////////////
    {
      size_t maxsamples 
	= std::min(magnet::GL::detail::glGet<GL_MAX_COLOR_TEXTURE_SAMPLES>(),
		   magnet::GL::detail::glGet<GL_MAX_DEPTH_TEXTURE_SAMPLES>());
      _aasamples.reset(new Gtk::ComboBoxText);
      for (size_t samples = maxsamples; samples > 0; samples /= 2)
	_aasamples->insert(0, boost::lexical_cast<std::string>(samples));
      
      _aasamples->set_active(0);
      _aasamples->show();
      
      Gtk::Box* AAsamplebox;
      _refXml->get_widget("AAbox", AAsamplebox);
      AAsamplebox->pack_start(*_aasamples, false, false);      
      
      _aasamples->signal_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::AAsamplechangeCallback));
      
      {
	Gtk::Entry* AmbientLightIntensity;
	_refXml->get_widget("AmbientLightIntensity", AmbientLightIntensity);
	AmbientLightIntensity
	  ->set_text(boost::lexical_cast<std::string>(_ambientIntensity));
	AmbientLightIntensity->signal_changed()
	  .connect(sigc::mem_fun(*this, &CLGLWindow::guiUpdateCallback));
      }
    
      {///////////////////////Filters//////////////////////////////////
	///Tree view must be built
      
	//Build the store
	_filterStore = Gtk::ListStore::create(*_filterModelColumns);
      
	//Setup the filter store
	_refXml->get_widget("filterView", _filterView);
	_filterView->set_model(_filterStore);
	_filterView->append_column("Active", _filterModelColumns->m_active);
	_filterView->append_column("Filter Name", _filterModelColumns->m_name);
      
	//////Connect the filterView select callback
	{
	  Glib::RefPtr<Gtk::TreeSelection> treeSelection
	    = _filterView->get_selection();
	
	  treeSelection->signal_changed()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterSelectCallback));
	}
      
	{///Connect the control buttons
	  Gtk::ToggleButton* btn;
	  _refXml->get_widget("filterActive", btn);
	  btn->signal_toggled()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterActiveCallback));
	}
      
	{
	  Gtk::CheckButton* btn;
	  _refXml->get_widget("filterEnable", btn);
	  btn->signal_toggled()
	    .connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
	}
      
	{//Fill the selector widgit with the available filters
	  Gtk::ComboBox* filterSelectBox;
	  _refXml->get_widget("filterSelectBox", filterSelectBox);
	  Filter::populateComboBox(filterSelectBox);
	}
      }
    
      {/////////////////////3D effects
	{
	  Gtk::CheckButton* stereoEnable;
	  _refXml->get_widget("StereoModeEnable", stereoEnable);
	  stereoEnable->signal_toggled()
	    .connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
	}

	{
	  Gtk::ScrolledWindow* win;
	  _refXml->get_widget("OpenVRLogScrolledWindow", win);	  	  
	  Gtk::TextView*  view;
	  _refXml->get_widget("OpenVRLogTextView", view);
	  view->signal_size_allocate().connect([=](Gtk::Allocation&) {
	      win->get_vadjustment()->set_value(win->get_vadjustment()->get_upper());
	    });

	  auto vrlog = Glib::RefPtr<Gtk::TextBuffer>::cast_dynamic(_refXml->get_object("OpenVRTextBuffer"));
	  Gtk::Button*  clear;
	  _refXml->get_widget("OpenVRLogClearButton", clear);
	  clear->signal_clicked().connect([=](){ vrlog->set_text(""); });
	}	
#ifdef COIL_OpenVR
	{
	  Gtk::CheckButton* vrEnable;
	  _refXml->get_widget("OpenVREnable", vrEnable);
	  vrEnable->signal_toggled().connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
	  vrEnable->set_sensitive(true);

	  auto vrlog = Glib::RefPtr<Gtk::TextBuffer>::cast_dynamic(_refXml->get_object("OpenVRTextBuffer"));
	  vrlog->set_text("OpenVR support available.\n");
	}
#else
	{
	  auto vrlog = Glib::RefPtr<Gtk::TextBuffer>::cast_dynamic(_refXml->get_object("OpenVRTextBuffer"));
	  vrlog->set_text("OpenVR not available (was not compiled in).\n");
	}
#endif
	
	{
	  Gtk::ComboBox* stereoMode;
	  _refXml->get_widget("StereoMode", stereoMode);
	  stereoMode->set_active(0);
	}

	
	{
	  Gtk::Entry* simunits;
	  _refXml->get_widget("SimLengthUnits", simunits);

	  std::ostringstream os;
	  os << _camera.getRenderScale();
	  simunits->set_text(os.str());

	  simunits->signal_changed()
	    .connect(sigc::bind(&magnet::gtk::forceNumericEntry, simunits));
	  simunits->signal_activate()
	    .connect(sigc::mem_fun(*this, &CLGLWindow::guiUpdateCallback));
	}
	
	{
	  Gtk::Entry* pixelPitch;
	  _refXml->get_widget("pixelPitch", pixelPitch);

	  std::ostringstream os;
	  os << _camera.getPixelPitch() * 10;
	  pixelPitch->set_text(os.str());

	  pixelPitch->signal_changed()
	    .connect(sigc::bind(&magnet::gtk::forceNumericEntry, pixelPitch));
	  pixelPitch->signal_activate()
	    .connect(sigc::mem_fun(*this, &CLGLWindow::guiUpdateCallback));
	}

#ifdef COIL_wiimote
	{//Here all the wii stuff should go in
	  Gtk::Button* btn;
	  _refXml->get_widget("wiiConnectBtn", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::wiiMoteConnect));
	  btn->set_sensitive(true);
	}
	
	{
	  Gtk::DrawingArea *ir;
	  _refXml->get_widget("wiiIRImage", ir);
	  ir->signal_draw().connect(sigc::mem_fun(this, &CLGLWindow::wiiMoteIRExposeEvent));	  
	}
	
	{//Here all the wii stuff should go in
	  Gtk::Button* btn;
	  _refXml->get_widget("wiiCalibrate", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(&(magnet::TrackWiimote::getInstance()), 
				   &magnet::TrackWiimote::calibrate));
	}
#endif
      }
    }

    {///////////////////////Render Objects//////////////////////////////////
      ///Tree view must be built    
      //Setup the filter store
      { 
	Gtk::TreeView* tree;
	_refXml->get_widget("renderObjView", tree);
	_renderObjsTree.init(tree);
      }
      
      //Populate the render object view
      _renderObjsTree.buildRenderView();
      selectRObjCallback();

      //////Connect the view to the select callback
      {
	Glib::RefPtr<Gtk::TreeSelection> treeSelection
	  = _renderObjsTree._view->get_selection();
      
	treeSelection->signal_changed()
	  .connect(sigc::mem_fun(this, &CLGLWindow::selectRObjCallback));
      }
    }

    if (_dynamo)
      {
	{
	  Gtk::Box* dynamoOpts;
	  _refXml->get_widget("dynamoOpts", dynamoOpts);
	
	  dynamoOpts->set_visible();
	}
      
	{
	  Gtk::Label* dynamoLabel;
	  _refXml->get_widget("simOptionsLabel", dynamoLabel);
	
	  dynamoLabel->set_visible();
	}

	{
	  Gtk::CheckButton* btn;
	  _refXml->get_widget("forceParticleSync", btn);
	  btn->signal_toggled()
	    .connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));

	  _particleSync = btn->get_active();
	}
      }
  }

  bool 
  CLGLWindow::GTKTick()
  {
    //This callback is used to calculate the FPS and sim update rates
    int currFrameTime = glutGet(GLUT_ELAPSED_TIME);

    float fps = _frameCounter * 1000.0 / (currFrameTime - _FPStime);
    float ups = _updateCounter * 1000.0 / (currFrameTime - _FPStime);
  
    std::stringstream fpsstring;
    fpsstring << "FPS:" << fps;
  
    Gtk::Label* label;
    _refXml->get_widget("RenderUpdateLabel", label);
    label->set_text(fpsstring.str());

    fpsstring.str("");
    fpsstring << "UPS:" << ups;

    _refXml->get_widget("SimUpdateLabel", label);
    label->set_text(fpsstring.str());

    _frameCounter = 0;
    _updateCounter = 0;
    _FPStime = currFrameTime;

    return true;
  }

  void
  CLGLWindow::init()
  {
    std::lock_guard<std::mutex> lock(_destroyLock);

    if (_readyFlag) return;

    const Vector look_at = Vector{0, 0, 0};
    const Vector up = Vector{0,1,0};
    {
      std::shared_ptr<RLight> light(new RLight("Light", Vector{1, 1, 1}, look_at, 0.1f, 300.0f, up, _camera.getRenderScale(), 0.2));
      _renderObjsTree._renderObjects.push_back(light);
    }
  
    _consoleID = _renderObjsTree._renderObjects.size();
    std::array<GLfloat, 3> textcolor  = {{0.5, 0.5, 0.5}};
    std::shared_ptr<RenderObj> consoleObj(new Console(textcolor)); 
    _renderObjsTree._renderObjects.push_back(consoleObj);

    glutInitContextVersion(3, 2);
    glutInitContextProfile(GLUT_CORE_PROFILE);
#ifdef MAGNET_DEBUG
    glutInitContextFlags(GLUT_DEBUG);
#endif
    glutInitDisplayMode(GLUT_RGBA);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(0, 0);

    CoilRegister::getCoilInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

    _glContext = magnet::GL::Context::getContext();

    glDepthFunc(GL_LEQUAL);
    _glContext->setDepthTest(true);

    //Setup the viewport
    _camera.resize(800, 600, _samples);
    _cairo_screen.init(800, 600);
 
    //Setup the keyboard controls
    glutIgnoreKeyRepeat(1);

    _lastUpdateTime = _lastFrameTime = _FPStime = glutGet(GLUT_ELAPSED_TIME);
   
    _copyShader.build();
    _downsampleShader.build();
    _blurShader.build();
    _pointLightShader.build();
    _shadowLightShader.build();
    _ambientLightShader.build();
    _luminanceShader.build();
    _luminanceMipMapShader.build();
    _toneMapShader.build();
    _depthResolverShader.build();
    
    {
      //Build depth buffer
      std::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D());
      //We don't force GL_DEPTH_COMPONENT24 as it is likely you get
      //the best precision anyway
      depthTexture->init(1024, 1024, GL_DEPTH_COMPONENT);//SIZE MUST BE THE SAME FOR THE LIGHTS
      //You must select GL_NEAREST for depth data, as GL_LINEAR
      //converts the value to 8bit for interpolation (on NVidia).
      depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
      
      //Build color texture
      std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D());
      colorTexture->init(1024, 1024, GL_RG32F);//SIZE MUST BE THE SAME FOR THE LIGHTS
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

      _shadowbuffer.init();
      _shadowbuffer.attachTexture(colorTexture, 0);
      _shadowbuffer.attachTexture(depthTexture);
    }

      //Now init the render objects  
    for (auto& obj: _renderObjsTree._renderObjects)
      obj->init(_systemQueue);
  
    try {
      initGTK();
    } catch(std::exception& err)
      {
	M_throw() << "An exception was thrown while initialising GTK\n"
		  << err.what() << "\n";
      }
    
    _readyFlag = true;
  }

  void
  CLGLWindow::deinit()
  {
    std::lock_guard<std::mutex> lock(_destroyLock);
  
    if (!_readyFlag) return;
    _readyFlag = false;

    ////////////////////GTK
    //Get rid of any filters, if we call the callback, a dialog will be instanced
    for (auto& child : _filterStore->children())
      {
	void* tmp_ptr = child[_filterModelColumns->m_filter_ptr];
	delete static_cast<Filter*>(tmp_ptr);
      }
    _filterStore->clear();

    _timeout_connection.disconnect();

    {
      Gtk::Window* controlwindow;
      _refXml->get_widget("controlWindow", controlwindow);  
      controlwindow->hide();
    }
  
    _refXml.reset(); //Destroy GTK instance
    _aasamples.reset();

    ///////////////////OpenGL
    for (auto& obj : _renderObjsTree._renderObjects) obj->deinit();

    _renderObjsTree._renderObjects.clear();
    _camera.deinit();    
    _shadowbuffer.deinit();
    _toneMapShader.deinit();
    _depthResolverShader.deinit();
    _pointLightShader.deinit();	
    _shadowLightShader.deinit();	
    _ambientLightShader.deinit();
    _downsampleShader.deinit();
    _blurShader.deinit();
    _copyShader.deinit();
    _luminanceShader.deinit();
    _luminanceMipMapShader.deinit();

    _cairo_screen.deinit();
#ifdef MAGNET_FFMPEG_SUPPORT
    _encoder.reset();
#endif
    ///////////////////Finally, unregister with COIL
    CoilRegister::getCoilInstance().unregisterWindow(this);
  }

  void 
  CLGLWindow::CallBackDisplayFunc()
  {
    if (!CoilRegister::getCoilInstance().isRunning() || !_readyFlag) return;
    //Setup the timings
    int _currFrameTime = glutGet(GLUT_ELAPSED_TIME);

    ////////////CAMERA UPDATES
    {
      Gtk::ToggleButton* button;
      _refXml->get_widget("CamRescalebtn", button);
      if (button->get_active())
	rescaleCameraCallback();
    }

    _camera.setRotatePoint(_cameraFocus);
    if (_selectedObject && (_cameraMode == ROTATE_POINT))
      {
	std::array<GLfloat, 4> vec = _selectedObject->getCursorPosition(_selectedObjectID);
	_camera.setRotatePoint(magnet::math::Vector{vec[0], vec[1], vec[2]});
      }

    float moveAmp  = (_currFrameTime - _lastFrameTime) * _moveSensitivity;

    float forward  = moveAmp * (keyStates[static_cast<size_t>('w')] - keyStates[static_cast<size_t>('s')]);
    float sideways = moveAmp * (keyStates[static_cast<size_t>('d')] - keyStates[static_cast<size_t>('a')]);
    float vertical =  moveAmp * (keyStates[static_cast<size_t>('q')] - keyStates[static_cast<size_t>('z')]);
    _camera.movement(0, 0, forward, sideways, vertical);

#ifdef COIL_wiimote
    //Run an update if the wiiMote was connected
    if ((magnet::TrackWiimote::getInstance()).connected())
      {
	Gtk::CheckButton* wiiHeadTrack;
	_refXml->get_widget("wiiHeadTracking", wiiHeadTrack);
	
	if (wiiHeadTrack->get_active())
	  _camera.setEyeLocation(magnet::TrackWiimote::getInstance().getHeadPosition());
      }
#endif

    bool enable_2D = true;
#ifdef COIL_OpenVR
    if (_openVRMode) {
      enable_2D = false;
      _openVR.getPosesAndSync();

      _openVR.setEye(vr::Eye_Left);
      drawScene(_openVR, enable_2D);
      _openVR.submit();

      _openVR.setEye(vr::Eye_Right);
      drawScene(_openVR, enable_2D);
      _openVR.submit();

      _openVR.PostPresentHandoff();
      
      _openVR.handleEvents();
    }
#endif

    ////////All of the camera movement and orientation has been
    ////////calculated with a certain fixed head position, now we
    ////////actually perform the rendering with adjustments for the
    ////////eyes
    
    const Vector oldHeadPosition = _camera.getEyeLocation();
    Vector headPosition = oldHeadPosition;

    if (!_stereoMode) {
      _camera.setEyeLocation(headPosition);
      drawScene(_camera, enable_2D);
      _camera._renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight());
    } else {
      const double eyedist = 6.5;
      Vector eyeDisplacement{0.5 * eyedist, 0, 0};
      
	Gtk::ComboBox* stereoMode;
	_refXml->get_widget("StereoMode", stereoMode);
	int mode = stereoMode->get_active_row_number();

	switch(mode)
	  {
	  case 0: //Analygraph Red-Cyan
	    //Do the right eye
	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera, enable_2D);

	    glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);	
	    _glContext->setDepthTest(false);
	    _camera._renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    _copyShader["u_Texture0"] = 0;
	    _copyShader.invoke(); 
	    _copyShader.detach();

	    ////Do the left eye
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera, enable_2D);

	    glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE);
	    _glContext->setDepthTest(false);
	    _camera._renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    _copyShader["u_Texture0"] = 0;
	    _copyShader.invoke(); 
	    _copyShader.detach();

	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	    break;
	  case 1:
	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera, enable_2D);
	    _camera._renderTarget.blitToScreen(_camera.getWidth() / 2, 
				       _camera.getHeight(), 0, 0, GL_LINEAR);
	    
	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera, enable_2D);
	    _camera._renderTarget.blitToScreen(_camera.getWidth() / 2, _camera.getHeight(),
				       _camera.getWidth() / 2, 0, GL_LINEAR);	    
	    break;
	  case 2:
	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera, enable_2D);
	    _camera._renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight()  /2,
				       0, 0, GL_LINEAR);
	    
	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera, enable_2D);
	    _camera._renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight() / 2,
				       0, _camera.getHeight() / 2, GL_LINEAR);
	    break;
	  default:
	    M_throw() << "Unknown stereo render mode";
	  }
      }
    //Reset the eye position
    _camera.setEyeLocation(oldHeadPosition);

    getGLContext()->swapBuffers();

    //Check if we're recording and then check that if we're
    //framelocking, check that new data is available
    if (_snapshot || ((_record) && (!_simframelock || _newData)))
      {	
	std::vector<uint8_t> pixels;
	pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	//Read the pixels into our container
	_camera._renderTarget.getColorTexture()->writeto(pixels);
	
	//Chop off the alpha channel to save bandwidth
	for (size_t i(0); i < _camera.getWidth() * _camera.getHeight(); ++i)
	  for (size_t component(0); component < 3; ++component)
	    pixels[3*i + component] = pixels[4 * i + component];
	
	pixels.resize(_camera.getWidth() * _camera.getHeight() * 3);

	std::string path;
	{
	  Gtk::FileChooserButton* fileChooser;
	  _refXml->get_widget("snapshotDirectory", fileChooser);
	  path = fileChooser->get_filename();
	}
	
	std::ostringstream filename;
	filename << std::setw(6) <<  std::setfill('0') << std::right << std::dec << _snapshot_counter++;
	
	if (_snapshot)
	  magnet::image::writePNGFile(path + "/" + filename.str() +".png", pixels, 
				      _camera.getWidth(), _camera.getHeight(), 3, 1, true, true);

#ifdef MAGNET_FFMPEG_SUPPORT	
	if (_record) _encoder->addFrame(pixels, true);
#endif
	_snapshot = false;
	_newData = false;
      }

    ++_frameCounter;
    _lastFrameTime = _currFrameTime;

    ////////////GUI UPDATES
    //We frequently ping the gui update
    guiUpdateCallback();
  }

  void 
  CLGLWindow::drawScene(magnet::GL::Camera& camera, bool draw_2D_overlay)
  {
    magnet::GL::FBO& renderTarget = camera.getResolveBuffer();
    
    //We perform a deffered shading pass followed by a forward shading
    //pass for objects which cannot be deferred, like volumes etc.

    ///////////////////////Deferred Shading G-Buffer Creation /////////////////
    //We use the stencil buffer to track which pixels should be shaded
    //in the deferred pass.

    //We share the depth and stencil texture between the GBuffer and
    //the target fbo
    camera._Gbuffer.attach();
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    _glContext->setDepthTest(true);
    _glContext->setBlend(false);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Enter the render ticks for all objects
    for (auto& obj :_renderObjsTree._renderObjects)
      if (obj->visible()) obj->glRender(camera, RenderObj::DEFAULT);

    camera._Gbuffer.detach();
    
    ///////////////////////Lighting pass////////////////////////
    //Here we calculate the lighting of every pixel in the scene
    camera._Gbuffer.getColorTexture(0)->bind(0);
    camera._Gbuffer.getColorTexture(1)->bind(1);
    camera._Gbuffer.getColorTexture(2)->bind(2);

    //First, set up the buffers for rendering
    camera._hdrBuffer.attach();
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //We need the depth test on, to enable writes to the depth buffer
    _glContext->setDepthTest(true);

    //Now we need to populate the depth buffer, but nothing needs to
    //go in the color buffer
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    _depthResolverShader.attach();
    _depthResolverShader["posTex"] = 2;
    _depthResolverShader["samples"] = GLint(_samples);
    _depthResolverShader["ProjectionMatrix"] = camera.getProjectionMatrix();
    _depthResolverShader.invoke();
    _depthResolverShader.detach();
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

    //Additive blending of all of the lights contributions, except for
    //the alpha values
    _glContext->setBlend(true);
    glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ZERO);

    //Now disable writing or testing of the depth buffer
    _glContext->setDepthTest(false);
    glDepthMask(GL_FALSE);

    _ambientLightShader.attach();
    _ambientLightShader["colorTex"] = 0;
    _ambientLightShader["samples"] = GLint(_samples);
    _ambientLightShader["ambientLight"] = GLfloat(_ambientIntensity);
    _ambientLightShader.invoke();
    _ambientLightShader.detach();

    std::vector<std::shared_ptr<RLight> > lights;

    //Perform the shadow casting lights
    for (auto& light_obj :_renderObjsTree._renderObjects)
      {
	std::shared_ptr<RLight> light = std::dynamic_pointer_cast<RLight>(light_obj);
	if (!light || !(light->shadowCasting())) continue;

	//Change from the hdr FBO 
	camera._hdrBuffer.detach();
	//Render each light's shadow map
	_shadowbuffer.attach();
	_glContext->setDepthTest(true);
	glDepthMask(GL_TRUE);
	_glContext->setBlend(false);

	glClearColor(light->getZFar(), light->getZFar() *light->getZFar(), 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	for (auto& obj :_renderObjsTree._renderObjects)
	  if (obj->visible() && obj->shadowCasting())
	    obj->glRender(*light, RenderObj::SHADOW);

	_shadowbuffer.detach();

	//_shadowbuffer.getColorTexture(0)->genMipmaps();
	_shadowbuffer.getColorTexture(0)->bind(7);

	//Change back to the hdr FBO
	camera._hdrBuffer.attach();
	_glContext->setDepthTest(false);
	glDepthMask(GL_FALSE);
	_glContext->setBlend(true);
	_shadowLightShader.attach();
	_shadowLightShader["colorTex"] = 0;
	_shadowLightShader["normalTex"] = 1;
	_shadowLightShader["positionTex"] = 2;
	_shadowLightShader["shadowTex"] = 7;
	_shadowLightShader["shadowMatrix"]
	  = light->getShadowTextureMatrix() * inverse(camera.getViewMatrix());
	_shadowLightShader["samples"] = GLint(_samples);
	_shadowLightShader["lightColor"] = light->getLightColor();
	_shadowLightShader["lightSpecularExponent"] = light->getSpecularExponent();
	_shadowLightShader["lightSpecularFactor"] = light->getSpecularFactor();
	_shadowLightShader["lightPosition"] = light->getEyespacePosition(camera);
	_shadowLightShader["maxVariance"] = light->getMaxVariance();
	_shadowLightShader["bleedReduction"] = light->getBleedReduction();
	_shadowLightShader.invoke();
	_shadowLightShader.detach();
      }

    //Perform the point/non-shadowing lights
    _pointLightShader.attach();
    _pointLightShader["colorTex"] = 0;
    _pointLightShader["normalTex"] = 1;
    _pointLightShader["positionTex"] = 2;
    _pointLightShader["samples"] = GLint(_samples);
    for (auto& light_obj :_renderObjsTree._renderObjects)
      {
	std::shared_ptr<RLight> light = std::dynamic_pointer_cast<RLight>(light_obj);
	if (!light || light->shadowCasting()) continue;
	lights.push_back(light);
	_pointLightShader["lightColor"] = light->getLightColor();
	_pointLightShader["lightSpecularExponent"] = light->getSpecularExponent();
	_pointLightShader["lightSpecularFactor"] = light->getSpecularFactor();
	_pointLightShader["lightPosition"] = light->getEyespacePosition(camera);
	_pointLightShader.invoke();
      }
    _pointLightShader.detach();

    _glContext->setBlend(true);
    _glContext->setDepthTest(true);
    glDepthMask(GL_TRUE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //Enter the forward render ticks for all objects
    for (auto& obj : _renderObjsTree._renderObjects)
      if (obj->visible())
	obj->forwardRender(camera._hdrBuffer, camera, lights, 
			   _ambientIntensity, RenderObj::DEFAULT);
    
    camera._hdrBuffer.detach();	
    ///////////////////////Luminance Sampling//////////////////////
    //The light buffer is bound to texture unit 0 for the tone mapping too
    _glContext->setDepthTest(false);
    _glContext->setBlend(false);

    camera._hdrBuffer.getColorTexture()->bind(0);

    camera._luminanceBuffer1.attach();
    _luminanceShader.attach();
    _luminanceShader["colorTex"] = 0;
    _luminanceShader.invoke();
    _luminanceShader.detach();
    camera._luminanceBuffer1.detach();

//    std::vector<GLfloat> data;
//    _luminanceBuffer1.getColorTexture()->writeto(data);
//
//    double avg = 0, max = data[1], min = data[2], count = 0;
//    for (size_t i(0); i < data.size() / 4; ++i)
//      {
//	if (data[4*i+3])
//	  {
//	    if (!count)
//	      { max = data[4*i+1]; min = data[4*i+2]; }
//	    
//	    avg += data[4*i+0] * data[4*i+3];
//	    max = std::max(double(data[4*i+1]), max);
//	    min = std::min(double(data[4*i+2]), min);
//	    count += data[4*i+3];
//	  }
//      }
//
//    avg /= count;
//    std::cout << "Luminance Avg="<< avg << ", max=" << max << ", min=" << min << ", count=" << count * 10000.0 << "\n";
      

    magnet::GL::FBO* luminanceSource = &camera._luminanceBuffer1;
    magnet::GL::FBO* luminanceDestination = &camera._luminanceBuffer2;

    //Now we need to generate the mipmaps containing the scene
    //average, minimum and maximum luminances.
    {
      GLsizei currentWidth = camera._luminanceBuffer1.getColorTexture()->getWidth();
      GLsizei currentHeight = camera._luminanceBuffer1.getColorTexture()->getHeight();
      GLint numLevels = camera._luminanceBuffer1.getColorTexture()->calcMipmapLevels();

      //Attach the mipmapping shader
      _luminanceMipMapShader.attach();
      for (int i=1; i < numLevels; ++i)
	{
	  luminanceDestination->attach();
	  luminanceSource->getColorTexture()->bind(0);
	  _luminanceMipMapShader["inputTex"] = 0;

	  std::array<GLint,2> oldSize = {{currentWidth, currentHeight}};
	  _luminanceMipMapShader["oldSize"] = oldSize;

	  //Halve the size of the textures, ensuring they never drop below 1
	  currentWidth /= 2; currentWidth += !currentWidth;
	  currentHeight /= 2; currentHeight += !currentHeight;
	  _glContext->setViewport(0, 0, currentWidth, currentHeight);

	  //Now generate the mipmap level using a shader
	  _luminanceMipMapShader.invoke();

	  luminanceDestination->detach();
	  std::swap(luminanceSource, luminanceDestination);
	}
      _luminanceMipMapShader.detach();
    }

    ///////////////////////Blurred Scene///////////////////////////////
    if (_bloomEnable)
      {
	magnet::GL::Texture2D& tex = *camera._hdrBuffer.getColorTexture();
	tex.bind(0);
      
	camera._blurTarget1.attach();
	_downsampleShader.attach();
	_downsampleShader["inputTex"] = 0;
	_downsampleShader["downscale"] = GLint(4);
	std::array<GLint,2> oldSize = {{tex.getWidth(), tex.getHeight()}};
	_downsampleShader["oldSize"] = oldSize;
	_downsampleShader.invoke();
	_downsampleShader.detach();
	camera._blurTarget1.detach();


	_blurShader.attach();
	_blurShader["colorTex"] = 0;
	std::array<GLfloat, 2> invDim = {{1.0f / (tex.getWidth() / 4),
					       1.0f / (tex.getHeight() / 4)}};
	_blurShader["invDim"] = invDim;

	for (size_t passes(0); passes < 1; ++passes)
	  {
	    camera._blurTarget1.getColorTexture()->bind(0);
	    camera._blurTarget2.attach();
	    _blurShader["direction"] = 0;	 
	    _blurShader.invoke();
	    _blurShader.detach();	  
	    camera._blurTarget2.detach();
	  
	    camera._blurTarget2.getColorTexture()->bind(0);
	    camera._blurTarget1.attach();
	    _blurShader.attach();
	    _blurShader["direction"] = 1;
	    _blurShader.invoke();
	    camera._blurTarget1.detach();
	  }
	_blurShader.detach();
      }

    ///////////////////////Tone Mapping///////////////////////////
    renderTarget.attach();
    camera._hdrBuffer.getColorTexture()->bind(0);
    luminanceSource->getColorTexture()->bind(1);
    if (_bloomEnable)
      camera._blurTarget1.getColorTexture()->bind(2);
    _toneMapShader.attach();
    _toneMapShader["color_tex"] = 0;
    _toneMapShader["logLuma"] = 1;
    _toneMapShader["bloom_tex"] = 2;
    _toneMapShader["bloom_enable"] = _bloomEnable;
    _toneMapShader["bloomCompression"] = GLfloat(_bloomCompression);
    _toneMapShader["bloomCutoff"] = GLfloat(_bloomCutoff);
    _toneMapShader["Lpwhite"] = GLfloat(_bloomSaturation);
    _toneMapShader["scene_key"] = GLfloat(_sceneKey);
    _toneMapShader["background_color"] = _backColor;
    _toneMapShader.invoke();
    _toneMapShader.detach();
    renderTarget.detach();

    //////////////////////FILTERING////////////
    //Attempt to perform some filtering

    bool FBOalternate = false;
    magnet::GL::FBO* lastFBO = &renderTarget;
    
    if (_filterEnable)
      {
       	//Bind the original image to texture unit 0
       	renderTarget.getColorTexture(0)->bind(0);

	//We can attach the GBuffer textures, for the normals and the
	//positions.
	//
       	//Normals unit 1
       	camera._Gbuffer.getColorTexture(1)->bind(1);
       	//Screen space positions 2
       	camera._Gbuffer.getColorTexture(2)->bind(2);
         
       	for (auto& child : _filterStore->children())
       	  {
       	    void* filter_ptr = child[_filterModelColumns->m_filter_ptr];
       	    Filter& filter = *static_cast<Filter*>(filter_ptr);
       	  
       	    if (!(child[_filterModelColumns->m_active])) continue; //Only run active filters, skip to the next filter
       	    if (filter.type_id() == detail::filterEnum<FlushToOriginal>::val)
       	      {//Check if we're trying to flush the drawing
       		lastFBO->attach();
       		glActiveTextureARB(GL_TEXTURE0);
       		//Now copy the texture 
       		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, camera.getWidth(), camera.getHeight());
       		lastFBO->detach();
       	      }
       	    else
       	      {
       		lastFBO->getColorTexture()->bind(3);
       		//The last output goes into texture 3
       		if (FBOalternate)
       		  camera._filterTarget1.attach();
       		else
       		  camera._filterTarget2.attach();
       	      
       		filter.invoke(3, camera.getWidth(), camera.getHeight(), camera);
       	      
       		if (FBOalternate)
       		  { camera._filterTarget1.detach(); lastFBO = &camera._filterTarget1; }
       		else
       		  { camera._filterTarget2.detach(); lastFBO = &camera._filterTarget2; }
       	      
       		FBOalternate = !FBOalternate;
       	      }
       	  }
      }
       
    //////////////Interface draw////////////////////////
    //We need alpha blending for the overlays

    if (draw_2D_overlay) { 
      _glContext->setBlend(true);
      lastFBO->attach();
      //Enter the interface draw for all objects
      _cairo_screen.clear();
      
      _glContext->cleanupAttributeArrays();
      for (auto& obj : _renderObjsTree._renderObjects)
	obj->interfaceRender(camera, _cairo_screen);
      
      //Draw the cursor if an object is selected
      if (_selectedObject)
	{
	  std::array<GLfloat, 4> vec = _selectedObject->getCursorPosition(_selectedObjectID);
	  vec = camera.project(Vector{vec[0], vec[1], vec[2]});
	  _cairo_screen.drawCursor(vec[0], vec[1], 5);
	  _cairo_screen.drawTextBox(vec[0] + 5, vec[1] + 5, 
				    _selectedObject->getCursorText(_selectedObjectID), 
				    5);
	}

      _cairo_screen.syncCairoGL();
      _cairo_screen.glRender();
      lastFBO->detach();
      _glContext->setBlend(false);
    }
    

    //Check if we actually did something and copy the data to the
    //output FBO if needed
    if (lastFBO != &renderTarget)
      {
	lastFBO->attach();
	renderTarget.getColorTexture(0)->bind(0);
	glActiveTextureARB(GL_TEXTURE0);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, camera.getWidth(), 
			    camera.getHeight());
	lastFBO->detach();
      }

    _glContext->setDepthTest(true);
  }

  void CLGLWindow::CallBackReshapeFunc(int w, int h)
  {
    if (!CoilRegister::getCoilInstance().isRunning() || !_readyFlag) return;
    resizeRender(w, h);
  }

  void CLGLWindow::resizeRender(int w, int h)
  {
    //We cannot resize a window below a threshold
    if ((w < 4) || (h < 4)) return;

    if ((size_t(h) == _camera.getHeight()) && (size_t(w) == _camera.getWidth()))
      return; //Skip a null op

    _camera.resize(w, h, _samples);
    _cairo_screen.resize(w, h);
  }

  void 
  CLGLWindow::setWindowtitle(const std::string& newtitle) 
  { 
    windowTitle = newtitle;
    glutSetWindowTitle(windowTitle.c_str());
  }

  void 
  CLGLWindow::CallBackMouseFunc(int button, int state, int x, int y)
  {
    switch (button)
      {
      case GLUT_LEFT_BUTTON:

	if (state == GLUT_DOWN)
	  {
	    _oldMouseX = x;
	    _oldMouseY = y;
	  
	    _mouseState |= LEFTMOUSE;
	  }
	else
	  _mouseState &= ~LEFTMOUSE;
	break;
      case GLUT_RIGHT_BUTTON:
	if (state == GLUT_DOWN)
	  {
	    _oldMouseX = x;
	    _oldMouseY = y;
	  
	    _mouseState |= RIGHTMOUSE;

	    //Now perform a picking selection
	    performPicking(x,y);
	  }
	else
	  _mouseState &= ~RIGHTMOUSE;
	break;
      case GLUT_MIDDLE_BUTTON:
	if (state == GLUT_DOWN)
	  {
	    _oldMouseX = x;
	    _oldMouseY = y;
	  
	    _mouseState |= MIDDLEMOUSE;
	  }
	else
	  _mouseState &= ~MIDDLEMOUSE;
	break;
      case 3:
	if (state == GLUT_UP) _moveSensitivity *= 1.1;
	break;
      case 4:
	if (state == GLUT_UP) _moveSensitivity /= 1.1;
	break;
      default:
	break;
      }
  }

  void 
  CLGLWindow::CallBackMouseWheelFunc(int button, int dir, int x, int y)
  {
    if (dir > 0)
      _moveSensitivity *= 1.1;
    else
      _moveSensitivity /= 1.1;
  }

  void 
  CLGLWindow::CallBackMotionFunc(int x, int y)
  {
    float diffY = (y-_oldMouseY) * _mouseSensitivity;
    float diffX = (x-_oldMouseX) * _mouseSensitivity;

    if (_mouseState & LEFTMOUSE)
      _camera.movement(diffX, diffY, 0, 0, 0);

    if (_mouseState & RIGHTMOUSE)
      if (_selectedObject)
	{
	  //Get the direction of the mouse from the camera position.
	  magnet::math::Vector ray = _camera.unprojectToDirection(x, y);

	  //We're dragging a selected object (the picking occurs on
	  //right mouse button down). Calculate the current position
	  //of the cursor
	  std::array<GLfloat, 4> vec 
	    = _selectedObject->getCursorPosition(_selectedObjectID);
	  const magnet::math::Vector origin{vec[0], vec[1], vec[2]};
	  const magnet::math::Vector camdir = _camera.getCameraDirection();
	  const magnet::math::Vector campos = _camera.getPosition();
	  const magnet::math::Vector rij = origin - campos;
	  double obj_distance = (rij | camdir) / camdir.nrm();
	  
	  double ray_distance = (ray | camdir) / camdir.nrm();
	  
	  ray *= obj_distance / ray_distance;
	  const magnet::math::Vector cursor_pos = campos + ray;
	  _selectedObject->dragCallback(cursor_pos, _selectedObjectID);
	}
  
    _oldMouseX = x;
    _oldMouseY = y;
  }

  void 
  CLGLWindow::CallBackKeyboardFunc(unsigned char key, int x, int y)
  {
    keyStates[std::tolower(key)] = true;
    
    if (std::tolower(key) == 'f')
      glutFullScreenToggle();
  }

  void 
  CLGLWindow::CallBackKeyboardUpFunc(unsigned char key, int x, int y)
  {
    keyStates[std::tolower(key)] = false;
  }

  bool
  CLGLWindow::simupdateTick(double t)
  {
    if (!isReady()) return false;
    
    //A loop for framelocked rendering, this holds the simulation
    //until the last data update has been rendered.
    while (_simframelock && (_lastUpdateTime == getLastFrameTime()))
      {
	//Jump out without an update if the window has been killed
	if (!isReady()) return false;
	_systemQueue->drainQueue();
	
	//1ms delay to lower CPU usage while blocking, but not to
	//affect framelocked render rates
	timespec sleeptime;
	sleeptime.tv_sec = 0;
	sleeptime.tv_nsec = 1000000;
	nanosleep(&sleeptime, NULL);
      }

    //Update the simulation data.  Only update if the previous data
    //set has been rendered, or we're about to pause the simulation
    if ((_lastUpdateTime != getLastFrameTime()) || !_simrun)
      {
	std::lock_guard<std::mutex> lock(_destroyLock);
	if (!isReady()) return false;
	_updateDataSignal();
	_newData = true;
	
	std::ostringstream os;
	os << "t:" << t;        
	setSimStatus1(os.str());
	_lastUpdateTime = getLastFrameTime();
      }

    ++_updateCounter;

    //A loop for paused running, to hold the system at the current
    //frame.
    while (!_simrun)
      {
	//Jump out without an update if the window has been killed
	if (!isReady()) return false;
	_systemQueue->drainQueue();
	
	//1ms delay to lower CPU usage while blocking
	timespec sleeptime;
	sleeptime.tv_sec = 0;
	sleeptime.tv_nsec = 1000000;
	nanosleep(&sleeptime, NULL);
      }

    return true;
  }

  void 
  CLGLWindow::runCallback()
  { 
    Gtk::ToggleButton* togButton;
    _refXml->get_widget("SimRunButton", togButton);

    Gtk::Image* togButtonImage;
    _refXml->get_widget("SimRunButtonImage", togButtonImage);
  
    Gtk::StockID origimage;
    Gtk::IconSize origsize;
    togButtonImage->get_stock(origimage, origsize);

    //Set the icon depending on the state
    if ((_simrun = togButton->get_active()))
      togButtonImage->set(Gtk::StockID("gtk-media-pause"), origsize);
    else
      togButtonImage->set(Gtk::StockID("gtk-media-play"), origsize);
  }

  void 
  CLGLWindow::simFramelockControlCallback()
  {
    Gtk::ToggleButton* framelockButton;
    _refXml->get_widget("SimLockButton", framelockButton);

    _simframelock = framelockButton->get_active();
  }

  void 
  CLGLWindow::snapshotCallback()
  {
    _snapshot = true;
  }

  void 
  CLGLWindow::recordCallback()
  {
    Gtk::ToggleButton* recordButton;
    _refXml->get_widget("SimRecordButton", recordButton);

#ifdef MAGNET_FFMPEG_SUPPORT
    if (_encoder.get() == NULL) 
      _encoder.reset(new magnet::image::VideoEncoderFFMPEG);
    

    if (_record != recordButton->get_active())
      {//The record button has been toggled
	if (!_record)
	  {
	    std::string path;
	    {
	      Gtk::FileChooserButton* fileChooser;
	      _refXml->get_widget("snapshotDirectory", fileChooser);
	      path = fileChooser->get_filename();
	    }
	
	    std::ostringstream counterstr;
	    counterstr << std::setw(6) <<  std::setfill('0') << std::right << std::dec << _video_counter++;
	    std::string filename =  path + "/" + counterstr.str() +".mpg";
	    
	    _encoder->open(filename, _camera.getWidth(), _camera.getHeight());
	  }
	else
	  _encoder->close();
      }
#endif

    _record = recordButton->get_active();  

  }

  void 
  CLGLWindow::filterUpCallback()
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _filterView->get_selection();

    Gtk::TreeModel::iterator iter_1 = refTreeSelection->get_selected();  
    Gtk::TreeModel::iterator iter_2 = iter_1;
    --iter_2;
    _filterStore->iter_swap(iter_1, iter_2);

    filterSelectCallback();
  }

  void 
  CLGLWindow::filterDownCallback()
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _filterView->get_selection();

    Gtk::TreeModel::iterator iter_1 = refTreeSelection->get_selected();  
    Gtk::TreeModel::iterator iter_2 = iter_1;
    ++iter_2;
    _filterStore->iter_swap(iter_1, iter_2);
  
    filterSelectCallback();
  }

  void 
  CLGLWindow::filterDeleteCallback()
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _filterView->get_selection();
  
    Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();
  
    void* tmp_ptr = (*iter)[_filterModelColumns->m_filter_ptr];
    delete static_cast<Filter*>(tmp_ptr);

    _filterStore->erase(iter);

    filterSelectCallback();
  }

  void 
  CLGLWindow::filterAddCallback()
  {
    //Grab the filter select box
    Gtk::ComboBox* filterSelectBox;
    _refXml->get_widget("filterSelectBox", filterSelectBox);

    //Check the filterSelectBox is on a valid row
    if (filterSelectBox->get_active_row_number() < 0) return;

    Gtk::TreeModel::iterator iter = _filterStore->append();

    size_t type_id = (*filterSelectBox->get_active())
      [Filter::getSelectColumnsInstance().m_col_id];

    (*iter)[_filterModelColumns->m_filter_ptr] 
      = Filter::createFromID(type_id);

    (*iter)[_filterModelColumns->m_name]
      = Filter::getName(type_id);
  
    (*iter)[_filterModelColumns->m_active]
      = true;

    filterSelectCallback();
  }

  void 
  CLGLWindow::filterSelectCallback()
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _filterView->get_selection();

    Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

    Gtk::Button *upbtn, *downbtn, *deletebtn;
    Gtk::ToggleButton *activeBtn;
    Gtk::Image *activeImage;
    _refXml->get_widget("filterUp", upbtn);
    _refXml->get_widget("filterDown", downbtn);
    _refXml->get_widget("filterDelete", deletebtn);
    _refXml->get_widget("filterActive", activeBtn);
    _refXml->get_widget("filterActiveImage", activeImage);

    Gtk::ScrolledWindow* frame;
    _refXml->get_widget("FilterOptions", frame);
    frame->remove();

    if(iter)
      {
	Gtk::TreeModel::iterator next_iter = iter;
	++next_iter;

	Filter* filter_ptr
	  = (Filter*)((void*)((*iter)[_filterModelColumns->m_filter_ptr]));
      
	//Enable the filter buttons
	upbtn    ->set_sensitive(iter != _filterStore->children().begin());
	downbtn  ->set_sensitive(next_iter);
	deletebtn->set_sensitive(true);
	activeBtn->set_sensitive(true);
      
	if (filter_ptr->getActive())
	  {//Object is visible
	    activeBtn->set_active(true);
	    activeImage->set(Gtk::Stock::YES, Gtk::ICON_SIZE_BUTTON);
	  }
	else
	  {//Object is not visible
	    activeBtn->set_active(false);
	    activeImage->set(Gtk::Stock::NO, Gtk::ICON_SIZE_BUTTON);
	  }

	filter_ptr->showControls(frame);
      }
    else
      {
	//Disable all of the filter buttons
	upbtn    ->set_sensitive(false);
	downbtn  ->set_sensitive(false); 
	deletebtn->set_sensitive(false);
	activeBtn ->set_sensitive(false);
      }
  }

  void
  CLGLWindow::filterActiveCallback()
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _filterView->get_selection();
    Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

    if (iter)
      {
	Gtk::ToggleButton *filterActive;
	_refXml->get_widget("filterActive", filterActive);
      
	bool newState = filterActive->get_active();
      
	Filter* filter_ptr
	  = (Filter*)((void*)((*iter)[_filterModelColumns->m_filter_ptr]));
	filter_ptr->setActive(newState);
	(*iter)[_filterModelColumns->m_active] = newState;
      }
  }

  void 
  CLGLWindow::SaveDataCallback() {
    Gtk::FileChooserDialog dialog(*controlwindow, "Save as", Gtk::FILE_CHOOSER_ACTION_SAVE);    
    Glib::RefPtr<Gtk::FileFilter> filter = Gtk::FileFilter::create();
    filter->set_name("Coil file");
    filter->add_pattern("*.coil");
    dialog.add_filter(filter);
    dialog.add_button("Ok", Gtk::RESPONSE_OK);
    dialog.add_button("Cancel", Gtk::RESPONSE_CANCEL);
    if (dialog.run() == Gtk::RESPONSE_OK) {
      std::string file_path = dialog.get_filename();
      if ((file_path.size() < 5) || (file_path.substr(file_path.size() - 5, 5) != ".coil"))
	file_path = file_path + ".coil";
      
      stator::xml::Document doc;
      auto root = doc.add_node("Coil");
      for (auto& obj :_renderObjsTree._renderObjects)
	obj->xml(root);
      std::ofstream(file_path, std::ios::out | std::ios::trunc) << doc;
    }
  }
  
  void 
  CLGLWindow::LoadDataCallback() {
    Gtk::FileChooserDialog dialog(*controlwindow, "Load data");

    std::vector<std::pair<std::string, std::string>> types{
      {"Coil file", "coil"},
      {"Volume raw", "raw"}
    };

    Glib::RefPtr<Gtk::FileFilter> filter = Gtk::FileFilter::create();
    filter->set_name("All supported types");
    for (const auto& pair : types)
      filter->add_pattern("*."+pair.second);
    dialog.add_filter(filter);
    
    for (const auto& pair : types) {
      Glib::RefPtr<Gtk::FileFilter> filter = Gtk::FileFilter::create();
      filter->set_name(pair.first + " (*."+pair.second+")");
      filter->add_pattern("*."+pair.second);
      dialog.add_filter(filter);
    }

    dialog.add_button("Ok", Gtk::RESPONSE_OK);
    dialog.add_button("Cancel", Gtk::RESPONSE_CANCEL);

    if (dialog.run() == Gtk::RESPONSE_OK) {
      std::string path = dialog.get_filename();
      std::string filename_only = boost::filesystem::path(dialog.get_filename()).filename().string();

      if ((filename_only.size() >= 4) && (filename_only.substr(filename_only.size()-4,4) == ".raw")) {
	Gtk::Dialog* volumedialog;
	dialog.hide();
	_refXml->get_widget("VolumeLoadDialog", volumedialog);
	volumedialog->set_transient_for(*controlwindow);
	
	std::ifstream in(dialog.get_filename(), std::ifstream::ate | std::ifstream::binary);
	size_t file_size = in.tellg();
	in.close();      
	
	Gtk::Label* label;
	_refXml->get_widget("VolumeFileSizeLabel", label);
	label->set_text(std::to_string(file_size));
	
	_refXml->get_widget("VolumeFileNameLabel", label);
	label->set_text(filename_only);
	
	Gtk::SpinButton* but;
	_refXml->get_widget("VolumeDataSizeButton", but);
	
	
	if (volumedialog->run() == Gtk::RESPONSE_OK){
	  std::shared_ptr<coil::RVolume> voldata(new coil::RVolume(filename_only));
	  addRenderObj(voldata);
	  
	  size_t data_size;
	  std::array<size_t, 3> data_dims;
	  _refXml->get_widget("VolumeDataSizeButton", but);
	  data_size = size_t(but->get_value_as_int());
	  
	  _refXml->get_widget("VolumeXDataSizeButton", but);
	  data_dims[0] = size_t(but->get_value_as_int());
	  _refXml->get_widget("VolumeYDataSizeButton", but);
	  data_dims[1] = size_t(but->get_value_as_int());
	  _refXml->get_widget("VolumeZDataSizeButton", but);
	  data_dims[2] = size_t(but->get_value_as_int());
	  
	  getGLContext()->queueTask(std::bind(&coil::RVolume::loadRawFile, voldata.get(), dialog.get_filename(), data_dims, data_size, Vector{0,0,0}));
	}
	volumedialog->hide();
      } else if ((filename_only.size() >= 5) && (filename_only.substr(filename_only.size()-5,5) == ".coil")) {
	
      } else {
	M_throw() << "Unhandled file type " << filename_only.substr(filename_only.size()-3,3);
      }
    }
  }
  
  void 
  CLGLWindow::filterClearCallback()
  {
    if (_filterStore->children().empty()) return;

    Gtk::Window* window;
    _refXml->get_widget("controlWindow", window);

    Gtk::MessageDialog confirmation(*window, "Are you sure you wish to erase all filters?",
				    false, Gtk::MESSAGE_QUESTION, Gtk::BUTTONS_OK_CANCEL, true);

    switch(confirmation.run())
      {
      case Gtk::RESPONSE_OK:
	{
	  for (auto& child : _filterStore->children())
	    {
	      void* tmp_ptr = child[_filterModelColumns->m_filter_ptr];
	      delete static_cast<Filter*>(tmp_ptr);
	    }
	
	  _filterStore->clear();
	}
      case Gtk::RESPONSE_CANCEL:
	break;
      default:
	M_throw() << "Unexpected return value!";
      }
  }

  void
  CLGLWindow::FPSLimitCallback()
  {
    Gtk::ToggleButton* fpslockButton;
    _refXml->get_widget("FPSLimit", fpslockButton);
    _fpsLimit = fpslockButton->get_active();

    Gtk::SpinButton* fpsButton;
    _refXml->get_widget("FPSLimitVal", fpsButton);
    _fpsLimitValue = fpsButton->get_value();

    //_renderTimeout.disconnect();
    //if (!_fpsLimit)
    //  _renderTimeout = Glib::signal_idle().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc));
    //else if (_fpsLimitValue != 0)
    //  _renderTimeout = Glib::signal_timeout().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc), 1000 / _fpsLimitValue, Glib::PRIORITY_DEFAULT_IDLE);
  }

  void
  CLGLWindow::aboutCallback()
  {
    Gtk::Window* aboutWindow;
    _refXml->get_widget("aboutSplashWindow", aboutWindow);
    aboutWindow->show();
  }


  namespace {
    struct IterFinder
    {
      IterFinder(std::shared_ptr<RenderObj> selected, Gtk::TreeModelColumn<RenderObj*> col):
	_selected(selected),
	_col(col)
      {}

      bool check(const Gtk::TreeModel::iterator& iter)
      {
	RenderObj* obj = (*iter)[_col];

	if (obj==_selected.get())
	  {
	    _iter = iter;
	    return true;
	  }

	return false;
      }

      Gtk::TreeModel::iterator _iter;
      std::shared_ptr<RenderObj> _selected;
      Gtk::TreeModelColumn<RenderObj*> _col;
    };
  }

  void
  CLGLWindow::performPicking(int x, int y)
  {
    _camera._renderTarget.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Perform unique coloring of screen objects, note that the value 0 is no object picked
    uint32_t offset = 1;
    //Now render the scene
    //Enter the render ticks for all objects
    for (auto& obj : _renderObjsTree._renderObjects)
      {
	const uint32_t n_objects = obj->pickableObjectCount();
	
	//If there are pickable objects and they are visible, then render them.
	if (n_objects)
	  {
	    obj->glRender(_camera, RenderObj::PICKING, offset);
	    offset += n_objects;
	  }
      }

    unsigned char pixel[4];
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);  
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel);    
    _camera._renderTarget.detach();
    
    //For debugging the picking render
    //_renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight());
    //getGLContext()->swapBuffers();

    //Now let the objects know what was picked
    _selectedObject.reset();
    size_t _selectedObjectGlobalID = pixel[0] 
      + 256 * (pixel[1] + 256 * (pixel[2] + 256 * pixel[3]));

    offset = 1;
    for (auto& obj : _renderObjsTree._renderObjects)
      { 
	const uint32_t n_objects = obj->pickableObjectCount();
	
	if ((_selectedObjectGlobalID >= offset) && (_selectedObjectGlobalID - offset) < n_objects)
	  {
	    _selectedObjectID = _selectedObjectGlobalID - offset;
	    _selectedObject = obj->getPickedObject(_selectedObjectID, obj);
	    break;
	  }
	offset += n_objects;
      }

    if (_selectedObject)
      {
	IterFinder finder(_selectedObject, _renderObjsTree._columns->m_obj);

	_renderObjsTree._store->foreach_iter(sigc::mem_fun(&finder, &IterFinder::check));

	if (finder._iter)
	  _renderObjsTree._view->get_selection()->select(finder._iter);
      }
    else
      if (_cameraMode == ROTATE_POINT)
	_cameraMode = ROTATE_CAMERA;
  }

  void CLGLWindow::selectRObjCallback() 
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _renderObjsTree._view->get_selection();

    Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

    Gtk::ScrolledWindow* win;
    _refXml->get_widget("ObjectOptions", win);

    win->remove(); //Clear the current object controls
    if (iter)
      {
	//Load the controls for the window
	RenderObj* obj = (*iter)[_renderObjsTree._columns->m_obj];
	obj->showControls(win);
      }
  }

  void
  CLGLWindow::setUpdateRateUnitToSteps(size_t defaultsteps)
  {
    {//Sim Update Frequency Control
      Gtk::SpinButton* updateButton;
      _refXml->get_widget("updateFreq", updateButton);
      updateButton->set_range(1,100000);
      updateButton->set_digits(0);
      updateButton->set_value(defaultsteps);
    }
  }

  void
  CLGLWindow::cameraModeCallback()
  {
    switch (_cameraMode)
      {
      case ROTATE_CAMERA:
	_cameraMode = ROTATE_WORLD;
	_camera.setMode(magnet::GL::CameraHeadTracking::ROTATE_POINT);
	break;
      case ROTATE_WORLD:
	if (_selectedObject)
	  {
	    _cameraMode = ROTATE_POINT;
	    _camera.setMode(magnet::GL::CameraHeadTracking::ROTATE_POINT);
	    break;
	  }
      case ROTATE_POINT:
	_cameraMode = ROTATE_CAMERA;
	_camera.setMode(magnet::GL::CameraHeadTracking::ROTATE_CAMERA);
	break;
      default:
	M_throw() << "Cannot change camera mode as it's in an unknown mode";
      }
  }

  void
  CLGLWindow::addLightCallback()
  {
    std::shared_ptr<RLight> light(new RLight("Light", Vector{0, 1, 0} * 50 / _camera.getRenderScale(), Vector{0, 0, 0}, 0.1f, 300.0f, Vector{0,1,0}, _camera.getRenderScale()));
    _renderObjsTree._renderObjects.push_back(light);
    _renderObjsTree._renderObjects.back()->init(_systemQueue);
    _renderObjsTree.buildRenderView();
  }

  void
  CLGLWindow::addFunctionCallback()
  {
    std::shared_ptr<RSurface> function(new RSurface("Function"));
    _renderObjsTree._renderObjects.push_back(function);
    _renderObjsTree._renderObjects.back()->init(_systemQueue);
    _renderObjsTree.buildRenderView();
  }

  void
  CLGLWindow::guiUpdateCallback()
  {
    {//Dynamo particle sync checkbox
      Gtk::CheckButton* btn;
      _refXml->get_widget("forceParticleSync", btn);
    
      _particleSync = btn->get_active();
    }

    {
      Gtk::Image* icon;
      _refXml->get_widget("CamModeimg", icon);

      switch (_cameraMode)
	{
	case ROTATE_CAMERA:
	  icon->set(coil::images::cammode_fps());
	  break;
	case ROTATE_WORLD:
	  icon->set(coil::images::cammode_rotate_world());
	  break;
	case ROTATE_POINT:
	  icon->set(coil::images::cammode_rotate_cursor());
	  break;
	default:
	  M_throw() << "Cannot find the appropriate icon for the camera mode";
	}
    }

    {
      Gtk::ColorButton* colorButton;
      _refXml->get_widget("BackColor", colorButton);
      
      Gdk::Color color = colorButton->get_color();
      
      _backColor[0] = GLfloat(color.get_red()) / G_MAXUSHORT;
      _backColor[1] = GLfloat(color.get_green()) / G_MAXUSHORT;
      _backColor[2] = GLfloat(color.get_blue()) / G_MAXUSHORT;
    }    

    {
      Gtk::Entry* entry;
      _refXml->get_widget("SceneKeyEntry", entry);
      magnet::gtk::forceNumericEntry(entry);
      try {
	_sceneKey = boost::lexical_cast<double>(entry->get_text());
      } catch(...) {}
    }

    {
      Gtk::Entry* entry;
      _refXml->get_widget("BloomSaturationEntry", entry);
      magnet::gtk::forceNumericEntry(entry);
      try {
	_bloomSaturation = boost::lexical_cast<double>(entry->get_text());
      } catch(...) {}
    }

    {
      Gtk::Entry* entry;
      _refXml->get_widget("BloomCutoffEntry", entry);
      magnet::gtk::forceNumericEntry(entry);
      try {
	_bloomCutoff = boost::lexical_cast<double>(entry->get_text());
      } catch(...) {}

      Gtk::CheckButton* btn;
      _refXml->get_widget("bloomEnable", btn);
      _bloomEnable = btn->get_active();
    }

    {
      Gtk::Entry* entry;
      _refXml->get_widget("BloomCompressionEntry", entry);
      magnet::gtk::forceNumericEntry(entry);
      try {
	_bloomCompression = boost::lexical_cast<double>(entry->get_text());
      } catch(...) {}
    }

    {
      Gtk::Entry* ambientIntensityEntry;
      _refXml->get_widget("AmbientLightIntensity", ambientIntensityEntry);
      magnet::gtk::forceNumericEntry(ambientIntensityEntry);
      try {
	_ambientIntensity
	  = boost::lexical_cast<double>(ambientIntensityEntry->get_text());
      } catch(...) {}
    }

    {//Filter enable/disable
      Gtk::CheckButton* btn;
      _refXml->get_widget("filterEnable", btn);
    
      _filterEnable = btn->get_active();
    }

    {//Sim Update Frequency Control
      Gtk::Entry* updateFreq;
      _refXml->get_widget("updateFreq", updateFreq);
    
      magnet::gtk::forceNumericEntry(updateFreq);
      if (updateFreq->get_text()[0] == '-')
	{
	  std::string value(updateFreq->get_text());
	  value.erase(value.begin());
	  updateFreq->set_text(value);
	}

      try { _updateIntervalValue 
	  = boost::lexical_cast<float>(updateFreq->get_text()); }
      catch(...) {}
    }

    {//Analygraph work
      Gtk::CheckButton* btn;
      _refXml->get_widget("StereoModeEnable", btn);    
      _stereoMode = btn->get_active();
    }

#ifdef COIL_OpenVR
    {//OpenVR
      auto vrlog = Glib::RefPtr<Gtk::TextBuffer>::cast_dynamic(_refXml->get_object("OpenVRTextBuffer"));
      
      Gtk::CheckButton* btn;
      _refXml->get_widget("OpenVREnable", btn);
      const bool newMode = btn->get_active();

      if (newMode != _openVRMode) {
	if (newMode) {	  
	  //Init OpenVR
	  _openVR.setLog([=](std::string line){ vrlog->insert_at_cursor(line+"\n"); });
	  _openVR.init();
	  if (_openVR.initialised())
	    _openVRMode = true;
	  else
	    btn->set_active(false);
	} else {
	  _openVR.shutdown();
	  _openVRMode = false;
	  btn->set_active(false);
	}
      }
    }
#endif
    
    {
      Gtk::Entry* simunits;
      _refXml->get_widget("SimLengthUnits", simunits);
      std::string val = simunits->get_text();
      double setval = 0.0;
      try {
	setval = boost::lexical_cast<double>(val);
      } catch (...) {}

      if (setval <= 0) setval = 25;

      _camera.setRenderScale(setval);
    }

    {
      Gtk::Entry* pixelPitch;
      _refXml->get_widget("pixelPitch", pixelPitch);
      std::string val = pixelPitch->get_text();

      double setval = 0.0;
      try {
	setval = boost::lexical_cast<double>(val);
      } catch (...) {}

      if (setval <= 0) setval = 0.25;

      _camera.setPixelPitch(setval / 10);
    }

    {
      Gtk::Label* XHead;
      _refXml->get_widget("XHead", XHead);
      Gtk::Label* YHead;
      _refXml->get_widget("YHead", YHead);
      Gtk::Label* ZHead;
      _refXml->get_widget("ZHead", ZHead);
      std::ostringstream os;
      os << _camera.getEyeLocation()[0] << "cm";
      XHead->set_text(os.str());
      os.str("");
      os << _camera.getEyeLocation()[1] << "cm";
      YHead->set_text(os.str());
      os.str("");
      os << _camera.getEyeLocation()[2] << "cm";
      ZHead->set_text(os.str());
    }

#ifdef COIL_wiimote
    {  
      Gtk::Label* statuslabel;
      _refXml->get_widget("wiiStatus", statuslabel);

      Gtk::Label* anglelabel;
      _refXml->get_widget("wiiAngleStatus", anglelabel);

      Gtk::ProgressBar* batteryBar;
      _refXml->get_widget("wiiBattery", batteryBar);

      Gtk::Button* wiiCalibrate;
      _refXml->get_widget("wiiCalibrate", wiiCalibrate);

      Gtk::DrawingArea *ir;
      _refXml->get_widget("wiiIRImage", ir);

      Gtk::Label* wiiXHead;
      _refXml->get_widget("wiiXHead", wiiXHead);
      Gtk::Label* wiiYHead;
      _refXml->get_widget("wiiYHead", wiiYHead);
      Gtk::Label* wiiZHead;
      _refXml->get_widget("wiiZHead", wiiZHead);

      Gtk::CheckButton* wiiHeadTrack;
      _refXml->get_widget("wiiHeadTracking", wiiHeadTrack);

      if ((magnet::TrackWiimote::getInstance()).connected())
	{
	  statuslabel->set_text("WiiMote Connected");

	  std::ostringstream os;
	  os << (magnet::TrackWiimote::getInstance()).getCalibrationAngle();
	  anglelabel->set_text(os.str());

	  Vector headPos = (magnet::TrackWiimote::getInstance()).getHeadPosition();

	  os.str("");
	  os << headPos[0] << "cm";
	  wiiXHead->set_text(os.str());
	  os.str("");
	  os << headPos[1] << "cm";
	  wiiYHead->set_text(os.str());
	  os.str("");
	  os << headPos[2] << "cm";
	  wiiZHead->set_text(os.str());

	  batteryBar->set_fraction((magnet::TrackWiimote::getInstance()).getBatteryLevel());

	  wiiCalibrate->set_sensitive(true);
	  wiiHeadTrack->set_sensitive(true);
	  {
	    Glib::RefPtr<Gdk::Window> win = ir->get_window();
	    if (win)
	      {
		Gdk::Rectangle r(0, 0, ir->get_allocation().get_width(),
				 ir->get_allocation().get_height());
		win->invalidate_rect(r, false);
	      }
	  }
	}
      else
	{
	  statuslabel->set_text("WiiMote Disconnected");
	  anglelabel->set_text("N/A");
	  wiiXHead->set_text("-");
	  wiiYHead->set_text("-");
	  wiiZHead->set_text("-");
	  batteryBar->set_fraction(0);
	  wiiCalibrate->set_sensitive(false);
	  wiiHeadTrack->set_sensitive(false);
	}
    }
#endif

    magnet::math::Vector pos = _camera.getPosition();
    std::ostringstream os;
    os << "Coil visualizer (" << _camera.getWidth() << "," << _camera.getHeight() << "), Camera Pos (" << pos[0] << "," << pos[1] << "," << pos[2] << ")";
    setWindowtitle(os.str());
  }


  void 
  CLGLWindow::setSimStatus1(std::string status)
  {
    Gtk::Label* label;
    _refXml->get_widget("SimDataLabel1", label);

    CoilRegister::getCoilInstance().getTaskQueue()
      .queueTask(std::bind(&CLGLWindow::setLabelText, this, label, status));
  }

  void 
  CLGLWindow::setSimStatus2(std::string status)
  {
    Gtk::Label* label;
    _refXml->get_widget("SimDataLabel2", label);
  
    CoilRegister::getCoilInstance().getTaskQueue()
      .queueTask(std::bind(&CLGLWindow::setLabelText, this, label, status));
  }

  void 
  CLGLWindow::setLabelText(Gtk::Label* label, std::string text)
  {
    label->set_text(text);
  }

  void 
  CLGLWindow::wiiMoteConnect()
  {
#ifdef COIL_wiimote
    if ((magnet::TrackWiimote::getInstance()).connected())
      {
	guiUpdateCallback();
	return;
      }

    Gtk::Window* window;
    _refXml->get_widget("controlWindow", window);
    Gtk::MessageDialog confirmation(*window, "Place the WiiMote in discovery mode (hit the <b>1</b> &amp; <b>2</b> buttons together)\nThen hit Ok.",
				    true, Gtk::MESSAGE_INFO, Gtk::BUTTONS_OK, true);

    confirmation.run();
    (magnet::TrackWiimote::getInstance()).connect();
#endif
  }

  bool 
  CLGLWindow::wiiMoteIRExposeEvent(const Cairo::RefPtr<Cairo::Context>& cr)
  {
#ifdef COIL_wiimote
    Gtk::DrawingArea *ir;
    _refXml->get_widget("wiiIRImage", ir);

    cr->set_source_rgb(0, 0, 0);
    cr->set_line_width(1);
    
    //Draw the tracked sources with a red dot, but only if there are just two sources!
    
    size_t trackeddrawn = 2;
    for (const auto& irdata : magnet::TrackWiimote::getInstance().getSortedIRData())
      {
	cr->save();
	if (trackeddrawn-- > 0)
	  cr->set_source_rgb(1, 0, 0);
	
	float x = ir->get_allocation().get_width() * (1 - float(irdata.x) / CWIID_IR_X_MAX);
	float y = ir->get_allocation().get_height() * (1 - float(irdata.y) / CWIID_IR_Y_MAX) ;
	
	cr->translate(x, y);
	cr->arc(0, 0, irdata.size + 1, 0, 2 * M_PI);
	cr->fill();	    
	cr->restore();
      }
#endif
    return true;
  }

  void 
  CLGLWindow::AAsamplechangeCallback()
  {
    if (_aasamples.get() != 0)
      _samples = boost::lexical_cast<size_t>(_aasamples->get_active_text());

    _camera.resize(_camera.getWidth(), _camera.getHeight(), _samples);
  }
  

  void 
  CLGLWindow::autoscaleView()
  {
    _glContext->queueTask(std::bind(&CLGLWindow::rescaleCameraCallback, this));
  }

  void 
  CLGLWindow::rescaleCameraCallback()
  {
    magnet::math::Vector min{HUGE_VAL,HUGE_VAL,HUGE_VAL};
    magnet::math::Vector max{-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

    for (auto& obj : _renderObjsTree._renderObjects)
      {
	magnet::math::Vector child_max = obj->getMaxCoord();
	magnet::math::Vector child_min = obj->getMinCoord();
	
	for (size_t i(0); i < 3; ++i) {
	  min[i] = std::min(min[i], child_min[i]);
	  max[i] = std::max(max[i], child_max[i]);
	}
      }
    //Catch the exceptional case where there is nothing rendered

    if (std::isinf(min[0]) || std::isinf(min[1]) || std::isinf(min[2])
	|| std::isinf(max[0]) || std::isinf(max[1]) || std::isinf(max[2]))
      return;

    double maxdim = std::max(max[0] - min[0], std::max(max[1] - min[1], max[2] - min[2]));

    double oldScale = _camera.getRenderScale();
    double newScale = 10.0 / maxdim;
    magnet::math::Vector centre = 0.5 * (min + max); 
    magnet::math::Vector shift = centre - _cameraFocus;
	
    //Try to reset the camera, in-case its dissappeared to nan or inf.
    _camera.setPosition(_camera.getPosition() * oldScale / newScale + shift);
    _camera.setRenderScale(newScale);
	
    _cameraFocus = centre;
    _cameraMode = ROTATE_WORLD;
    _camera.setMode(magnet::GL::CameraHeadTracking::ROTATE_POINT);
    {
      Gtk::Entry* simunits;
      _refXml->get_widget("SimLengthUnits", simunits);
	  
      std::ostringstream os;
      os << _camera.getRenderScale();
      simunits->set_text(os.str());
    }
	
    //Shift the lighting for the scene
    for (auto& obj : _renderObjsTree._renderObjects)
      {
	std::shared_ptr<RLight> ptr 
	  = std::dynamic_pointer_cast<RLight>(obj);
	if (ptr)
	  {
	    ptr->setSize(ptr->getSize() * oldScale / newScale);
	    ptr->setIntensity(ptr->getIntensity() * (oldScale * oldScale) / (newScale * newScale));
	    ptr->setPosition(ptr->getPosition() * oldScale / newScale + shift);
	  }
      }
  }

  void 
  CLGLWindow::HeadReset()
  {
    //Assume the user is around 70cm from the screen
    _camera.setEyeLocation(Vector{0, 0, 70});
  }
}

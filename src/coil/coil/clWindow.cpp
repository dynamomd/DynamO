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
#include <magnet/function/task.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <gtkmm/volumebutton.h>
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
    _fpsLimit(true),
    _fpsLimitValue(25),
    _filterEnable(true),
    _stereoMode(false),
    _ambientIntensity(0.0005),
    _snapshot_counter(0),
    _video_counter(0),
    _samples(1),
    _dynamo(dynamo),
    _cameraMode(ROTATE_WORLD)
  {
    for (size_t i(0); i < 3; ++i)
      _backColor[i] = 1.0;

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
    _timeout_connection
      = Glib::signal_timeout().connect_seconds(sigc::mem_fun(this, &CLGLWindow::GTKTick), 1);

    //Timeout for render
    _renderTimeout = Glib::signal_timeout().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc), 
						    1000 / _fpsLimitValue, Glib::PRIORITY_DEFAULT_IDLE);

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
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(1,0,0)));

      _refXml->get_widget("CamPlusYbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(0,1,0)));

      _refXml->get_widget("CamPlusZbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(0,0,1)));


      _refXml->get_widget("CamNegXbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(-1,0,0)));

      _refXml->get_widget("CamNegYbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(0,-1,0)));

      _refXml->get_widget("CamNegZbtn", button);
      button->signal_clicked()
	.connect(sigc::bind(sigc::mem_fun(_camera, &magnet::GL::Camera::setViewAxis), magnet::math::Vector(0,0,-1)));
    }

    {
      Gtk::Button* button;    
      _refXml->get_widget("CamMode", button); 
      
      button->signal_clicked()
	.connect(sigc::mem_fun(*this, &CLGLWindow::cameraModeCallback));
    }

    {
      Gtk::Button* button;    
      _refXml->get_widget("addLightButton", button); 
      
      button->signal_clicked()
	.connect(sigc::mem_fun(*this, &CLGLWindow::addLightCallback));
    }
    
    {
      Gtk::Button* button;    
      _refXml->get_widget("addFunctionButton", button); 
      
      button->signal_clicked()
	.connect(sigc::mem_fun(*this, &CLGLWindow::addFunctionCallback));
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

    {//////Snapshot button
      Gtk::Button* btn;
      _refXml->get_widget("SimSnapshot", btn);
      btn->signal_clicked().connect(sigc::mem_fun(this, &CLGLWindow::snapshotCallback));    
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
	_aasamples->prepend(boost::lexical_cast<std::string>(samples));
      
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
	  Gtk::Button* btn;
	  _refXml->get_widget("filterUp", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterUpCallback));
	  _refXml->get_widget("filterDown", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterDownCallback));
	  _refXml->get_widget("filterDelete", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterDeleteCallback));
	  _refXml->get_widget("filterAdd", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterAddCallback));
	  _refXml->get_widget("filterClear", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::filterClearCallback));
	  {
	    Gtk::ToggleButton* btn;
	    _refXml->get_widget("filterActive", btn);
	    btn->signal_toggled()
	      .connect(sigc::mem_fun(this, &CLGLWindow::filterActiveCallback));
	  }    
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

	{
	  Gtk::Button* btn;
	  _refXml->get_widget("HeadTrackReset", btn);
	  btn->signal_clicked()
	    .connect(sigc::mem_fun(this, &CLGLWindow::HeadReset));
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
	  ir->signal_expose_event()
	    .connect(sigc::mem_fun(this, &CLGLWindow::wiiMoteIRExposeEvent));	  
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
    magnet::thread::ScopedLock lock(_destroyLock);
    
    if (_readyFlag) return;

    double light_distance = 50 / _camera.getRenderScale();
    Vector look_at = Vector(0.0f, 0.0f, 0.0f);
    Vector up = Vector(0,1,0);
    
    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 1", Vector(1, -1, 0) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 2", Vector(0, -1, 1) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 3", Vector(-std::sqrt(0.5), -1, -std::sqrt(0.5)) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 4", -Vector(1, -1, 0) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 5", -Vector(0, -1, 1) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 6", -Vector(-std::sqrt(0.5), -1, -std::sqrt(0.5)) * light_distance, look_at, 30.0, 10000.0f, up, _camera.getRenderScale()));
      _renderObjsTree._renderObjects.push_back(light);
    }
  
    _consoleID = _renderObjsTree._renderObjects.size();
    std::tr1::array<GLfloat, 3> textcolor  = {{0.5, 0.5, 0.5}};
    std::tr1::shared_ptr<RenderObj> consoleObj(new Console(textcolor)); 
    _renderObjsTree._renderObjects.push_back(consoleObj);

    glutInitContextVersion(3, 3);
    glutInitContextProfile(GLUT_CORE_PROFILE);
#ifdef MAGNET_DEBUG
    glutInitContextFlags(GLUT_DEBUG);
#endif
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(0, 0);

    CoilRegister::getCoilInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

    _glContext = magnet::GL::Context::getContext();

    glDepthFunc(GL_LEQUAL);
    _glContext->setDepthTest(true);

    //Setup the viewport
    resizeRender(800, 600);
 
    //Setup the keyboard controls
    glutIgnoreKeyRepeat(1);

    _lastUpdateTime = _lastFrameTime = _FPStime = glutGet(GLUT_ELAPSED_TIME);
    _frameRenderTime = 0;

    _copyShader.build();
    _downsampleShader.build();
    _blurShader.build();
    _pointLightShader.build();
    _ambientLightShader.build();
    _VSMShader.build();
    _luminanceShader.build();
    _luminanceMipMapShader.build();
    _toneMapShader.build();
    _depthResolverShader.build();
    
    _cairo_screen.init(600, 600);

      //Now init the render objects  
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->init(_systemQueue);
  
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
    magnet::thread::ScopedLock lock(_destroyLock);
  
    if (!_readyFlag) return;
    _readyFlag = false;

    ////////////////////GTK
    //Get rid of any filters, if we call the callback, a dialog will be instanced
    for (Gtk::TreeModel::iterator iPtr = _filterStore->children().begin();
	 iPtr; ++iPtr)
      {
	void* tmp_ptr = (*iPtr)[_filterModelColumns->m_filter_ptr];
	delete static_cast<Filter*>(tmp_ptr);
      }
    _filterStore->clear();

    _timeout_connection.disconnect();
    _renderTimeout.disconnect();

    {
      Gtk::Window* controlwindow;
      _refXml->get_widget("controlWindow", controlwindow);  
      controlwindow->hide_all();
    }
  
    _refXml.reset(); //Destroy GTK instance
    _aasamples.reset();

    ///////////////////OpenGL
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->deinit();

    _renderObjsTree._renderObjects.clear();

    _renderTarget.deinit();
    _Gbuffer.deinit();
    _hdrBuffer.deinit();
    _luminanceBuffer1.deinit();
    _luminanceBuffer2.deinit();
    _filterTarget1.deinit();
    _filterTarget2.deinit();
    _blurTarget1.deinit();
    _blurTarget2.deinit();
    _toneMapShader.deinit();
    _depthResolverShader.deinit();
    _pointLightShader.deinit();	
    _ambientLightShader.deinit();
    _VSMShader.deinit();
    _downsampleShader.deinit();
    _blurShader.deinit();
    _copyShader.deinit();
    _luminanceShader.deinit();
    _luminanceMipMapShader.deinit();

    _cairo_screen.deinit();
    _encoder.reset();
    ///////////////////Finally, unregister with COIL
    CoilRegister::getCoilInstance().unregisterWindow(this);
  }

  void 
  CLGLWindow::CallBackDisplayFunc()
  {
    if (!CoilRegister::getCoilInstance().isRunning()
	|| !_readyFlag) return;
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
	std::tr1::array<GLfloat, 4> vec = _selectedObject->getCursorPosition(_selectedObjectID);
	_camera.setRotatePoint(magnet::math::Vector(vec[0], vec[1], vec[2]));
      }

    float moveAmp  = (_currFrameTime - _lastFrameTime) * _moveSensitivity;

    float forward  = moveAmp * (keyStates[static_cast<size_t>('w')] - keyStates[static_cast<size_t>('s')]);
    float sideways = moveAmp * (keyStates[static_cast<size_t>('d')] - keyStates[static_cast<size_t>('a')]);
    float vertical =  moveAmp * (keyStates[static_cast<size_t>('q')] - keyStates[static_cast<size_t>('z')]);
    _camera.movement(0, 0, forward, sideways, vertical);

    ////////////GUI UPDATES
    //We frequently ping the gui update
    guiUpdateCallback();

    ////////////Lighting shadow map creation////////////////////
    //This stage only needs to be performed once per frame
    
//    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
//	   = _renderObjsTree._renderObjects.begin();
//	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
//      if ((*iPtr)->shadowCasting() && std::tr1::dynamic_pointer_cast<RLight>(*iPtr))
//	{
//	  std::tr1::shared_ptr<RLight> light = std::tr1::static_pointer_cast<RLight>(*iPtr);
//	  if (light)
//	    {
//	      _VSMShader.attach();
//	      //Render each light's shadow map
//	      _VSMShader["ProjectionMatrix"] = light->getProjectionMatrix();
//	      _VSMShader["ViewMatrix"] = light->getViewMatrix();	  
//	      light->shadowFBO().attach();
//	      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//	      //Enter the render ticks for all objects
//	      for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator jPtr 
//		     = _renderObjsTree._renderObjects.begin();
//		   jPtr != _renderObjsTree._renderObjects.end(); ++jPtr)
//		if (iPtr != jPtr)
//		  if ((*jPtr)->shadowCasting() && (*jPtr)->visible())
//		    (*iPtr)->glRender(*light, RenderObj::SHADOW);
//	
//	      light->shadowFBO().detach();
//	      /////////////MIPMAPPED shadow maps don't seem to work
//	      //_light0.shadowTex()->genMipmaps();
//	      light->shadowTex()->bind(7);
//	
//	      _VSMShader.detach();
//	    }
//	}
    ////////All of the camera movement and orientation has been
    ////////calculated with a certain fixed head position, now we
    ////////actually perform the rendering with adjustments for the 
    
    const Vector oldHeadPosition = _camera.getEyeLocation();
    Vector headPosition = oldHeadPosition;

#ifdef COIL_wiimote
    //Run an update if the wiiMote was connected
    if ((magnet::TrackWiimote::getInstance()).connected())
      {
	Gtk::CheckButton* wiiHeadTrack;
	_refXml->get_widget("wiiHeadTracking", wiiHeadTrack);
	
	if (wiiHeadTrack->get_active())
	  headPosition = magnet::TrackWiimote::getInstance().getHeadPosition();
      }
#endif

    //Bind to the multisample buffer
    if (!_stereoMode)
      {
	_camera.setEyeLocation(headPosition);
	drawScene(_camera);
	_renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight());
	_camera.setEyeLocation(oldHeadPosition);
      }
    else
      {
	const double eyedist = 6.5;
	Vector eyeDisplacement(0.5 * eyedist, 0, 0);

	Gtk::ComboBox* stereoMode;
	_refXml->get_widget("StereoMode", stereoMode);
	int mode = stereoMode->get_active_row_number();

	switch(mode)
	  {
	  case 0: //Analygraph Red-Cyan
	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    _copyShader["u_Texture0"] = 0;
	    glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
	    _copyShader.invoke(); 
	    _copyShader.detach();
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE);
	    _copyShader.invoke(); 
	    _copyShader.detach();
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	    break;
	  case 1:
	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.blitToScreen(_camera.getWidth() / 2, 
				       _camera.getHeight(), 0, 0, GL_LINEAR);

	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.blitToScreen(_camera.getWidth() / 2, _camera.getHeight(),
				       _camera.getWidth() / 2, 0, GL_LINEAR);	    
	    break;
	  case 2:
	    _camera.setEyeLocation(headPosition + eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight()  /2,
				       0, 0, GL_LINEAR);

	    _camera.setEyeLocation(headPosition - eyeDisplacement);
	    drawScene(_camera);
	    _renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight() / 2,
				       0, _camera.getHeight() / 2, GL_LINEAR);
	    break;
	  default:
	    M_throw() << "Unknown stereo render mode";
	  }
	//Reset the eye position
	_camera.setEyeLocation(oldHeadPosition);
      }

    getGLContext()->swapBuffers();
    glFinish();

    //Check if we're recording and then check that if we're
    //framelocking, check that new data is available
    if (_snapshot 
	|| ((_record) && (!_simframelock || _newData)))
      {	
	std::vector<uint8_t> pixels;
	pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	//Read the pixels into our container
	_renderTarget.getColorTexture()->writeto(pixels);
	
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
    _frameRenderTime = glutGet(GLUT_ELAPSED_TIME) - _currFrameTime;
  }

  void 
  CLGLWindow::drawScene(magnet::GL::Camera& camera)
  {
    //We perform a deffered shading pass followed by a forward shading
    //pass for objects which cannot be deferred, like volumes etc.

    ///////////////////////Deferred Shading G-Buffer Creation /////////////////
    //We use the stencil buffer to track which pixels should be shaded
    //in the deferred pass.

    //We share the depth and stencil texture between the GBuffer and
    //the target fbo
    _Gbuffer.attach();
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    _glContext->setDepthTest(true);
    _glContext->setBlend(false);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible()) 
	(*iPtr)->glRender(camera, RenderObj::DEFAULT);

    _Gbuffer.detach();
    
    ///////////////////////Lighting pass////////////////////////
    //Here we calculate the lighting of every pixel in the scene
    _Gbuffer.getColorTexture(0)->bind(0);
    _Gbuffer.getColorTexture(1)->bind(1);
    _Gbuffer.getColorTexture(2)->bind(2);

    //First, set up the buffers for rendering
    _hdrBuffer.attach();
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
    _depthResolverShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
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

    _pointLightShader.attach();
    _pointLightShader["colorTex"] = 0;
    _pointLightShader["normalTex"] = 1;
    _pointLightShader["positionTex"] = 2;
    _pointLightShader["samples"] = GLint(_samples);

    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if (std::tr1::dynamic_pointer_cast<RLight>(*iPtr))
	{
	  std::tr1::shared_ptr<RLight> light 
	    = std::tr1::dynamic_pointer_cast<RLight>(*iPtr);

	  _pointLightShader["lightColor"] = light->getLightColor();
	  _pointLightShader["lightSpecularExponent"] = light->getSpecularExponent();
	  _pointLightShader["lightSpecularFactor"] = light->getSpecularFactor();
	  _pointLightShader["lightPosition"] = light->getEyespacePosition(camera);
	  _pointLightShader.invoke();
	}
    
    _pointLightShader.detach();

    ///////////////////////Forward Shading Pass /////////////////
    std::vector<std::tr1::shared_ptr<RLight> > lights;
    

    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if (std::tr1::dynamic_pointer_cast<RLight>(*iPtr))
	lights.push_back(std::tr1::dynamic_pointer_cast<RLight>(*iPtr));

    _glContext->setBlend(true);
    _glContext->setDepthTest(true);
    glDepthMask(GL_TRUE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //Enter the forward render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible())
	(*iPtr)->forwardRender(_hdrBuffer, camera, lights, 
			       _ambientIntensity, RenderObj::DEFAULT);
    
    _hdrBuffer.detach();	
    ///////////////////////Luminance Sampling//////////////////////
    //The light buffer is bound to texture unit 0 for the tone mapping too
    _glContext->setDepthTest(false);
    _glContext->setBlend(false);

    _hdrBuffer.getColorTexture()->bind(0);

    _luminanceBuffer1.attach();
    _luminanceShader.attach();
    _luminanceShader["colorTex"] = 0;
    _luminanceShader.invoke();
    _luminanceShader.detach();
    _luminanceBuffer1.detach();

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
      

    magnet::GL::FBO* luminanceSource = &_luminanceBuffer1;
    magnet::GL::FBO* luminanceDestination = &_luminanceBuffer2;

    //Now we need to generate the mipmaps containing the scene
    //average, minimum and maximum luminances.
    {
      GLsizei currentWidth = _luminanceBuffer1.getColorTexture()->getWidth();
      GLsizei currentHeight = _luminanceBuffer1.getColorTexture()->getHeight();
      GLint numLevels = _luminanceBuffer1.getColorTexture()->calcMipmapLevels();

      //Attach the mipmapping shader
      _luminanceMipMapShader.attach();
      for (int i=1; i < numLevels; ++i)
	{
	  luminanceDestination->attach();
	  luminanceSource->getColorTexture()->bind(0);
	  _luminanceMipMapShader["inputTex"] = 0;

	  std::tr1::array<GLint,2> oldSize = {{currentWidth, currentHeight}};
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
	magnet::GL::Texture2D& tex = *_hdrBuffer.getColorTexture();
	tex.bind(0);
      
	_blurTarget1.attach();
	_downsampleShader.attach();
	_downsampleShader["inputTex"] = 0;
	_downsampleShader["downscale"] = GLint(4);
	std::tr1::array<GLint,2> oldSize = {{tex.getWidth(), tex.getHeight()}};
	_downsampleShader["oldSize"] = oldSize;
	_downsampleShader.invoke();
	_downsampleShader.detach();
	_blurTarget1.detach();


	_blurShader.attach();
	_blurShader["colorTex"] = 0;
	std::tr1::array<GLfloat, 2> invDim = {{1.0f / (tex.getWidth() / 4),
					       1.0f / (tex.getHeight() / 4)}};
	_blurShader["invDim"] = invDim;

	for (size_t passes(0); passes < 1; ++passes)
	  {
	    _blurTarget1.getColorTexture()->bind(0);
	    _blurTarget2.attach();
	    _blurShader["direction"] = 0;	 
	    _blurShader.invoke();
	    _blurShader.detach();	  
	    _blurTarget2.detach();
	  
	    _blurTarget2.getColorTexture()->bind(0);
	    _blurTarget1.attach();
	    _blurShader.attach();
	    _blurShader["direction"] = 1;
	    _blurShader.invoke();
	    _blurTarget1.detach();
	  }
	_blurShader.detach();
      }

    ///////////////////////Tone Mapping///////////////////////////
    _renderTarget.attach();
    _hdrBuffer.getColorTexture()->bind(0);
    luminanceSource->getColorTexture()->bind(1);
    if (_bloomEnable)
      _blurTarget1.getColorTexture()->bind(2);
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
    _renderTarget.detach();

    //////////////////////FILTERING////////////
    //Attempt to perform some filtering

    bool FBOalternate = false;
    magnet::GL::FBO* lastFBO = &_renderTarget;
    
    if (_filterEnable)
      {
       	//Bind the original image to texture unit 0
       	_renderTarget.getColorTexture(0)->bind(0);

	//We can attach the GBuffer textures, for the normals and the
	//positions.
	//
       	//Normals unit 1
       	_Gbuffer.getColorTexture(1)->bind(1);
       	//Screen space positions 2
       	_Gbuffer.getColorTexture(2)->bind(2);
         
       	for (Gtk::TreeModel::iterator iPtr = _filterStore->children().begin();
       	     iPtr != _filterStore->children().end(); ++iPtr)
       	  {
       	    void* filter_ptr = (*iPtr)[_filterModelColumns->m_filter_ptr];
       	    Filter& filter = *static_cast<Filter*>(filter_ptr);
       	  
       	    if (!((*iPtr)[_filterModelColumns->m_active])) continue; //Only run active filters, skip to the next filter
       	    if (filter.type_id() == detail::filterEnum<FlushToOriginal>::val)
       	      {//Check if we're trying to flush the drawing
       		lastFBO->attach();
       		glActiveTextureARB(GL_TEXTURE0);
       		//Now copy the texture 
       		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, _camera.getWidth(), _camera.getHeight());
       		lastFBO->detach();
       	      }
       	    else
       	      {
       		lastFBO->getColorTexture()->bind(3);
       		//The last output goes into texture 3
       		if (FBOalternate)
       		  _filterTarget1.attach();
       		else
       		  _filterTarget2.attach();
       	      
       		filter.invoke(3, _camera.getWidth(), _camera.getHeight(), _camera);
       	      
       		if (FBOalternate)
       		  { _filterTarget1.detach(); lastFBO = &_filterTarget1; }
       		else
       		  { _filterTarget2.detach(); lastFBO = &_filterTarget2; }
       	      
       		FBOalternate = !FBOalternate;
       	      }
       	  }
      }
       
    //////////////Interface draw////////////////////////
    //We need alpha blending for the overlays
    _glContext->setBlend(true);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    lastFBO->attach();
    //Enter the interface draw for all objects
    _cairo_screen.clear();

    _glContext->cleanupAttributeArrays();
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->interfaceRender(_camera, _cairo_screen);

    //Draw the cursor if an object is selected
    if (_selectedObject)
      {
	std::tr1::array<GLfloat, 4> vec = _selectedObject->getCursorPosition(_selectedObjectID);
	vec = camera.project(Vector(vec[0], vec[1], vec[2]));
	_cairo_screen.drawCursor(vec[0], vec[1], 5);
	_cairo_screen.drawTextBox(vec[0] + 5, vec[1] + 5, 
				  _selectedObject->getCursorText(_selectedObjectID), 
				  5);
      }

    _cairo_screen.syncCairoGL();
    _cairo_screen.glRender();
    lastFBO->detach();

    _glContext->setBlend(false);

    //Check if we actually did something and copy the data to the
    //output FBO if needed
    if (lastFBO != &_renderTarget)
      {
	lastFBO->attach();
	_renderTarget.getColorTexture(0)->bind(0);
	glActiveTextureARB(GL_TEXTURE0);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, _camera.getWidth(), 
			    _camera.getHeight());
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
    //We cannot resize a window 
    if ((w < 4) || (h < 4)) return;

    _camera.setHeightWidth(h, w);
    _renderTarget.deinit();
    _Gbuffer.deinit();
    _hdrBuffer.deinit();
    _luminanceBuffer1.deinit();
    _luminanceBuffer2.deinit();
    _filterTarget1.deinit();
    _filterTarget2.deinit();
    _blurTarget1.deinit();
    _blurTarget2.deinit();

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _filterTarget1.init();
      _filterTarget1.attachTexture(colorTexture, 0);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _filterTarget2.init();
      _filterTarget2.attachTexture(colorTexture, 0);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth() / 4, _camera.getHeight() / 4, GL_RGB16F);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _blurTarget1.init();
      _blurTarget1.attachTexture(colorTexture, 0);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth() / 4, _camera.getHeight() / 4, GL_RGB16F);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _blurTarget2.init();
      _blurTarget2.attachTexture(colorTexture, 0);
    }

    {
      //Build the main/left-eye render buffer
      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);

      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	depthTexture(new magnet::GL::Texture2D);
      depthTexture->init(_camera.getWidth(), _camera.getHeight(), 
			 GL_DEPTH_COMPONENT);
      depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

      _renderTarget.init();
      _renderTarget.attachTexture(colorTexture, 0);
      _renderTarget.attachTexture(depthTexture);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	depthTexture(new magnet::GL::Texture2D);
      depthTexture->init(_camera.getWidth(), _camera.getHeight(), 
			 GL_DEPTH_COMPONENT);
      depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

      _hdrBuffer.init();
      _hdrBuffer.attachTexture(colorTexture, 0);
      _hdrBuffer.attachTexture(depthTexture);
    }
      
    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	colorTexture(new magnet::GL::Texture2D);
	
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

      _luminanceBuffer1.init();
      _luminanceBuffer1.attachTexture(colorTexture, 0);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> 
	colorTexture(new magnet::GL::Texture2D);
	
      colorTexture->init(_camera.getWidth()/2, _camera.getHeight()/2, GL_RGBA16F);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

      _luminanceBuffer2.init();
      _luminanceBuffer2.attachTexture(colorTexture, 0);
    }
    AAsamplechangeCallback();

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
	  std::tr1::array<GLfloat, 4> vec 
	    = _selectedObject->getCursorPosition(_selectedObjectID);
	  const magnet::math::Vector origin(vec[0], vec[1], vec[2]);
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

  void
  CLGLWindow::simupdateTick(double t)
  {
    if (!isReady()) return;
    
    //A loop for framelocked rendering, this holds the simulation
    //until the last data update has been rendered.
    while (_simframelock && (_lastUpdateTime == getLastFrameTime()))
      {
	//Jump out without an update if the window has been killed
	if (!isReady()) return;
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
	magnet::thread::ScopedLock lock(_destroyLock);
	if (!isReady()) return;
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
	if (!isReady()) return;
	_systemQueue->drainQueue();
	
	//1ms delay to lower CPU usage while blocking
	timespec sleeptime;
	sleeptime.tv_sec = 0;
	sleeptime.tv_nsec = 1000000;
	nanosleep(&sleeptime, NULL);
      }
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
	  for (Gtk::TreeModel::iterator iPtr = _filterStore->children().begin();
	       iPtr; ++iPtr)
	    {
	      void* tmp_ptr = (*iPtr)[_filterModelColumns->m_filter_ptr];
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

    _renderTimeout.disconnect();
    if (!_fpsLimit)
      _renderTimeout = Glib::signal_timeout().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc), 10, 
						      Glib::PRIORITY_DEFAULT_IDLE);
    else if (_fpsLimitValue != 0)
      _renderTimeout = Glib::signal_timeout().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc), 
						      1000 / _fpsLimitValue, Glib::PRIORITY_DEFAULT_IDLE);
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
      IterFinder(std::tr1::shared_ptr<RenderObj> selected, Gtk::TreeModelColumn<RenderObj*> col):
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
      std::tr1::shared_ptr<RenderObj> _selected;
      Gtk::TreeModelColumn<RenderObj*> _col;
    };
  }

  void
  CLGLWindow::performPicking(int x, int y)
  {
    _renderTarget.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Perform unique coloring of screen objects, note that the value 0 is no object picked
    uint32_t offset = 1;
    //Now render the scene
    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator 
	   iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      {
	const uint32_t n_objects = (*iPtr)->pickableObjectCount();
	
	//If there are pickable objects and they are visible, then render them.
	if (n_objects)
	  {
	    (*iPtr)->pickingRender(_camera, offset);
	    offset += n_objects;
	  }
      }

    unsigned char pixel[4];
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);  
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel);    
    _renderTarget.detach();
    
    //For debugging the picking render
    //_renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight());
    //getGLContext()->swapBuffers();

    //Now let the objects know what was picked
    _selectedObject.reset();
    size_t _selectedObjectGlobalID = pixel[0] 
      + 256 * (pixel[1] + 256 * (pixel[2] + 256 * pixel[3]));

    offset = 1;
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator 
	   iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      { 
	const uint32_t n_objects = (*iPtr)->pickableObjectCount();
	
	if ((_selectedObjectGlobalID >= offset) && (_selectedObjectGlobalID - offset) < n_objects)
	  {
	    _selectedObjectID = _selectedObjectGlobalID - offset;
	    _selectedObject = (*iPtr)->getPickedObject(_selectedObjectID, *iPtr);
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
	_camera.setMode(magnet::GL::Camera::ROTATE_POINT);
	break;
      case ROTATE_WORLD:
	if (_selectedObject)
	  {
	    _cameraMode = ROTATE_POINT;
	    _camera.setMode(magnet::GL::Camera::ROTATE_POINT);
	    break;
	  }
      case ROTATE_POINT:
	_cameraMode = ROTATE_CAMERA;
	_camera.setMode(magnet::GL::Camera::ROTATE_CAMERA);
	break;
      default:
	M_throw() << "Cannot change camera mode as it's in an unknown mode";
      }
  }

  void
  CLGLWindow::addLightCallback()
  {
    std::tr1::shared_ptr<RLight> light(new RLight("Light",
						  Vector(0, 0, 0),
						  Vector(0, -1, 0),//Lookat
						  75.0f//Beam angle
						  ));
    _renderObjsTree._renderObjects.push_back(light);
    _renderObjsTree._renderObjects.back()->init(_systemQueue);
    _renderObjsTree.buildRenderView();
  }

  void
  CLGLWindow::addFunctionCallback()
  {
    std::tr1::shared_ptr<RSurface> function(new RSurface("Function"));
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
      .queueTask(magnet::function::Task::makeTask(&CLGLWindow::setLabelText, this, label, status));
  }

  void 
  CLGLWindow::setSimStatus2(std::string status)
  {
    Gtk::Label* label;
    _refXml->get_widget("SimDataLabel2", label);
  
    CoilRegister::getCoilInstance().getTaskQueue()
      .queueTask(magnet::function::Task::makeTask(&CLGLWindow::setLabelText, this, label, status));
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
  CLGLWindow::wiiMoteIRExposeEvent(GdkEventExpose* event)
  {
#ifdef COIL_wiimote
    Gtk::DrawingArea *ir;
    _refXml->get_widget("wiiIRImage", ir);

    Glib::RefPtr<Gdk::Window> window = ir->get_window();
    if (window)
      {
	Cairo::RefPtr<Cairo::Context> cr = window->create_cairo_context();

	if(event)
	  {
	    // clip to the area indicated by the expose event so that we only
	    // redraw the portion of the window that needs to be redrawn
	    cr->rectangle(event->area.x, event->area.y,
			  event->area.width, event->area.height);
	    cr->clip();
	  }
      
	cr->set_source_rgb(0, 0, 0);
	cr->set_line_width(1);

	//Draw the tracked sources with a red dot, but only if there are just two sources!
      
	const std::vector<magnet::TrackWiimote::IRData>& irdata 
	  = magnet::TrackWiimote::getInstance().getSortedIRData();
      
	size_t trackeddrawn = 2;
	for (std::vector<magnet::TrackWiimote::IRData>::const_iterator iPtr = irdata.begin();
	     iPtr != irdata.end(); ++iPtr)
	  {
	    cr->save();
	    if (trackeddrawn-- > 0)
	      cr->set_source_rgb(1, 0, 0);

	    float x = ir->get_allocation().get_width() * (1 - float(iPtr->x) / CWIID_IR_X_MAX);
	    float y = ir->get_allocation().get_height() * (1 - float(iPtr->y) / CWIID_IR_Y_MAX) ;

	    cr->translate(x, y);
	    cr->arc(0, 0, iPtr->size + 1, 0, 2 * M_PI);
	    cr->fill();	    
	    cr->restore();
	  }
      }
#endif
    return true;
  }

  void 
  CLGLWindow::AAsamplechangeCallback()
  {
    if (_aasamples.get() != 0)
      _samples = boost::lexical_cast<size_t>(_aasamples->get_active_text());

    //Build G buffer      
    std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2DMultisampled(_samples));
    colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
    
    std::tr1::shared_ptr<magnet::GL::Texture2D> normalTexture(new magnet::GL::Texture2DMultisampled(_samples));
    normalTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
    
    std::tr1::shared_ptr<magnet::GL::Texture2D> posTexture(new magnet::GL::Texture2DMultisampled(_samples));
    posTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
    
    std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2DMultisampled(_samples));
    depthTexture->init(_camera.getWidth(), _camera.getHeight(), GL_DEPTH_COMPONENT);    

    _Gbuffer.deinit();
    _Gbuffer.init();

    _Gbuffer.attachTexture(colorTexture, 0);
    _Gbuffer.attachTexture(normalTexture, 1);
    _Gbuffer.attachTexture(posTexture, 2);
    _Gbuffer.attachTexture(depthTexture);
  }
  

  void 
  CLGLWindow::autoscaleView()
  {
    _glContext->queueTask(magnet::function::Task::makeTask(&CLGLWindow::rescaleCameraCallback, this));
  }

  void 
  CLGLWindow::rescaleCameraCallback()
  {
    magnet::math::Vector min(HUGE_VAL,HUGE_VAL,HUGE_VAL);
    magnet::math::Vector max(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);

    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      {
	magnet::math::Vector child_max = (*iPtr)->getMaxCoord();
	magnet::math::Vector child_min = (*iPtr)->getMinCoord();
	
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
    double newScale = 40.0 / maxdim;
    magnet::math::Vector centre = 0.5 * (min + max); 
    magnet::math::Vector shift = centre - _cameraFocus;
	
    //Try to reset the camera, in-case its dissappeared to nan or inf.
    _camera.setPosition(_camera.getPosition() * oldScale / newScale + shift);
    _camera.setRenderScale(newScale);
	
    _cameraFocus = centre;
    _cameraMode = ROTATE_WORLD;
    _camera.setMode(magnet::GL::Camera::ROTATE_POINT);
    {
      Gtk::Entry* simunits;
      _refXml->get_widget("SimLengthUnits", simunits);
	  
      std::ostringstream os;
      os << _camera.getRenderScale();
      simunits->set_text(os.str());
    }
	
    //Shift the lighting for the scene
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      {
	std::tr1::shared_ptr<RLight> ptr 
	  = std::tr1::dynamic_pointer_cast<RLight>(*iPtr);
	if (ptr)
	  {
	    ptr->setSize(ptr->getSize() * oldScale / newScale);
	    ptr->setPosition(ptr->getPosition() * oldScale / newScale +shift);
	  }
      }
  }

  void 
  CLGLWindow::HeadReset()
  {
    //Assume the user is around 70cm from the screen
    _camera.setEyeLocation(Vector(0, 0, 70));
  }
}

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

#include <coil/clWindow.hpp>
#include <coil/RenderObj/Function.hpp>
#include <coil/RenderObj/console.hpp>
#include <coil/RenderObj/Volume.hpp>
#include <coil/RenderObj/Light.hpp>

#include <magnet/GL/context.hpp>

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

//The glade xml file is "linked" into a binary file and stuffed in the
//executable, these are the symbols to its data
extern const char _binary_clwingtk_gladexml_start[];
extern const char _binary_clwingtk_gladexml_end[];

//All of the icons and images compiled into the exectuable
extern const guint8 coilicon[];
extern const size_t coilicon_size;

extern const guint8 coilsplash[];
extern const size_t coilsplash_size;

extern const guint8 camplusx[];
extern const size_t camplusx_size;
extern const guint8 camplusy[];
extern const size_t camplusy_size;
extern const guint8 camplusz[];
extern const size_t camplusz_size;
extern const guint8 camnegx[];
extern const size_t camnegx_size;
extern const guint8 camnegy[];
extern const size_t camnegy_size;
extern const guint8 camnegz[];
extern const size_t camnegz_size;

extern const guint8 cammode_rotate[];
extern const size_t cammode_rotate_size;

extern const guint8 cammode_fps[];
extern const size_t cammode_fps_size;

namespace coil {
  CLGLWindow::CLGLWindow(std::string title,
			 double updateIntervalValue,
			 bool dynamo
			 ):
    _systemQueue(new magnet::thread::TaskQueue),
    _updateIntervalValue(updateIntervalValue),
    keyState(DEFAULT),
    windowTitle(title),
    _frameCounter(0),
    _updateCounter(0),
    _mouseSensitivity(0.3),
    _moveSensitivity(0.001),
    _specialKeys(0),
    _simrun(false),
    _simframelock(false),
    _snapshot(false),
    _record(false),
    _PNGFileFormat(true),
    _fpsLimit(true),
    _fpsLimitValue(35),
    _filterEnable(true),
    _stereoMode(false),
    _burnoutFactor(0.8),
    _ambientIntensity(0.01),
    //This value of the exposure is actually the subjective midde
    //brightness (18% reflectance) of photographic paper
    _exposure(0.18),
    _snapshot_counter(0),
    _samples(1),
    _dynamo(dynamo)
  {
    for (size_t i(0); i < 3; ++i)
      _backColor[i] = 1.0;

    for (size_t i(0); i < 256; ++i) keyStates[i] = false;
  }

  CLGLWindow::~CLGLWindow() {}

  bool
  CLGLWindow::CallBackIdleFunc()
  {
    try {
      glutSetWindow(windowID);
      CallBackDisplayFunc();
    } catch (cl::Error err)
      {
	std::cerr << "\n Window render caught an OpenCL exception\n"
		  << "An OpenCL error occured," << err.what()
		  << "\nError num of " << err.err()
		  << "\n As we're in a thread we can only exit(1)!";
	std::exit(1);
      } catch (std::exception& except)
      {
	std::cerr << "\n Window render caught a std::exception\n"
		  << except.what();
	std::exit(1);
      }  catch (...)
      {
	std::cerr << "\nRender thread caught an unknown exception!\n";
	std::exit(1);
      }

    return true;
  }

  void
  CLGLWindow::initGTK()
  {
    _filterModelColumns.reset(new FilterModelColumnsType);

    try
      {
	_refXml = Gtk::Builder::create_from_string
	  (std::string(_binary_clwingtk_gladexml_start,
		       _binary_clwingtk_gladexml_end));  
      }
    catch (std::exception& err)
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
    controlwindow->set_icon(Gdk::Pixbuf::create_from_inline
			    (coilicon_size, coilicon));

    {
      Gtk::Image* icon;

      _refXml->get_widget("CamPlusXimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camplusx_size, camplusx));

      _refXml->get_widget("CamPlusYimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camplusy_size, camplusy));

      _refXml->get_widget("CamPlusZimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camplusz_size, camplusz));

      _refXml->get_widget("CamNegXimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camnegx_size, camnegx));

      _refXml->get_widget("CamNegYimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camnegy_size, camnegy));

      _refXml->get_widget("CamNegZimg", icon);
      icon->set(Gdk::Pixbuf::create_from_inline(camnegz_size, camnegz));
    }

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

    ///////Register the about button
    {
      Gtk::ImageMenuItem* aboutButton;
      _refXml->get_widget("aboutItem", aboutButton);

      aboutButton->signal_activate()
	.connect(sigc::mem_fun(this, &CLGLWindow::aboutCallback));
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
      AAsamplechangeCallback();
      
      _aasamples->signal_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::AAsamplechangeCallback));

      {
	Gtk::Entry* exposureEntry;
	_refXml->get_widget("ExposureEntry", exposureEntry);
	exposureEntry->set_text(boost::lexical_cast<std::string>(_exposure));
      }
      
      {
	Gtk::Entry* burnoutEntry;
	_refXml->get_widget("BurnoutEntry", burnoutEntry);
	burnoutEntry->signal_changed()
	  .connect(sigc::mem_fun(*this, &CLGLWindow::guiUpdateCallback));
      }

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
	  os << _camera.getSimUnitLength();
	  simunits->set_text(os.str());

	  simunits->signal_changed()
	    .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *simunits));
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
	    .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *pixelPitch));
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
    
    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 1",
							   Vector(-0.8f,  0.8f, 0.8f),//Position
							   Vector(0.0f, 0.0f, 0.0f),//Lookat
							   75.0f//Beam angle
							   ));
      _renderObjsTree._renderObjects.push_back(light);
    }

    {
      std::tr1::shared_ptr<RLight> light(new RLight("Light 2",
							   Vector(0.8f,  0.8f, -0.8f),//Position
							   Vector(0.0f, 0.0f, 0.0f),//Lookat
							   75.0f//Beam angle
							   ));
      _renderObjsTree._renderObjects.push_back(light);
    }
  
    //First render object is the ground
    std::tr1::shared_ptr<RenderObj> groundObj
      (new RFunction((size_t)64,
		     Vector(-5, -0.6, -5),
		     Vector(10,0,0), Vector(0,0,10), Vector(0,1,0), //Axis of the function, x,y,z
		     -1, -1,//Start point of the functions evaluation (x,y)
		     1, 1,//Range of the function to evaluate (xrange,yrange
		     false, //Render a set of Axis as well?
		     true, //Is the shape static, i.e. is there no time dependence
		     "Ground",
		     "f=0;\n",
		     "normal = (float4)(0,0,1,0);\n",
		     "colors[0] = (uchar4)(255,255,255,255);"
		     ));

    _renderObjsTree._renderObjects.push_back(groundObj);

    //Second render object is the console
    _consoleID = _renderObjsTree._renderObjects.size();
    std::tr1::array<GLfloat, 3> textcolor  = {{0.5, 0.5, 0.5}};
    std::tr1::shared_ptr<RenderObj> consoleObj(new Console(textcolor)); 
    _renderObjsTree._renderObjects.push_back(consoleObj);

    glutInitContextVersion(3, 3);
#ifdef MAGNET_DEBUG
    glutInitContextFlags(GLUT_CORE_PROFILE | GLUT_DEBUG);
#else
    glutInitContextFlags(GLUT_CORE_PROFILE);
#endif
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(0, 0);

    CoilRegister::getCoilInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

    _glContext = magnet::GL::Context::getContext();

    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);

    //Setup the viewport
    CallBackReshapeFunc(800, 600);

    //Setup the keyboard controls
    glutIgnoreKeyRepeat(1);

    _lastUpdateTime = _lastFrameTime = _FPStime = glutGet(GLUT_ELAPSED_TIME);
    _frameRenderTime = 0;

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _filterTarget1.init();
      _filterTarget1.attachTexture(colorTexture, 0);
    }

    {
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _filterTarget2.init();
      _filterTarget2.attachTexture(colorTexture, 0);
    }

    _renderShader.build();
    _copyShader.build();
    _pointLightShader.build();
    _VSMShader.build();
    _simpleRenderShader.build();
    _luminanceShader.build();
    _luminanceMipMapShader.build();
    _toneMapShader.build();

    {
      {
	//Build the main/left-eye render buffer
	std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
	colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	//We use a shared depth/stencil buffer for the deferred and forward shading passes
	std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
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

	colorTexture->init(_camera.getWidth()/2, _camera.getHeight()/2, GL_RGB16F);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	_lightBuffer.init();
	_lightBuffer.attachTexture(colorTexture, 0);
      }
      
      {
	std::tr1::shared_ptr<magnet::GL::Texture2D> 
	  colorTexture(new magnet::GL::Texture2D);
	
	colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RG16F);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	colorTexture->genMipmaps(); //Ensure the mipmap chain is built/available

	_luminanceBuffer.init();
	_luminanceBuffer.attachTexture(colorTexture, 0);
      }

      {
	//Build depth buffer
	std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
	depthTexture->init(1024, 1024, GL_DEPTH_COMPONENT);
	depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
	
	//Build color texture
	std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
	colorTexture->init(1024, 1024, GL_RG32F);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	
	_shadowBuffer.init();
	_shadowBuffer.attachTexture(colorTexture, 0);
	_shadowBuffer.attachTexture(depthTexture);
      }
    }



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

    /////////////////OpenCL

    getGLContext()->getCLCommandQueue().finish();

    ///////////////////OpenGL
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->deinit();

    _renderObjsTree._renderObjects.clear();

    _renderTarget.deinit();
    _Gbuffer.deinit();
    _lightBuffer.deinit();
    _luminanceBuffer.deinit();
    _shadowBuffer.deinit();
	
    _filterTarget1.deinit();
    _filterTarget2.deinit();
    _renderShader.deinit();
    _toneMapShader.deinit();
    _pointLightShader.deinit();	
    _VSMShader.deinit();
    _simpleRenderShader.deinit();
    _copyShader.deinit();
    _luminanceShader.deinit();
    _luminanceMipMapShader.deinit();
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

    //Prepare for the OpenCL ticks
    glFinish();
//    {
//      //Add the sync and flush it into the device
//      GLsync waitFence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
//      glFlush();
//
//      //Now wait on the sync
//      for (;;)
//	{
//	  switch (glClientWaitSync(waitFence, 0, 100))
//	    {
//	    case GL_TIMEOUT_EXPIRED:
//	      //We are still waiting, perform some Gtk tasks
//	      //Gtk::Main::iteration(false);
//	      continue; //test again
//	    case GL_WAIT_FAILED:
//	      M_throw() << "Failed to sync the OpenGL queue";
//	    case GL_CONDITION_SATISFIED:
//	    case GL_ALREADY_SIGNALED:
//	      //These cases indicate the sync has occurred
//	      break;
//	    }
//	  break;
//	}
//    }

    //Run every objects OpenCL stage
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->clTick(_camera);

    //Camera Positioning

    float moveAmp  = (_currFrameTime - _lastFrameTime) * _moveSensitivity;      
    float forward  = moveAmp * ( keyStates[static_cast<size_t>('w')] 
				 - keyStates[static_cast<size_t>('s')]);
    float sideways = moveAmp * ( keyStates[static_cast<size_t>('d')] 
				 - keyStates[static_cast<size_t>('a')]);
    float vertical = moveAmp * ( keyStates[static_cast<size_t>('q')] 
				 - keyStates[static_cast<size_t>('z')]);
    _camera.CameraUpdate(forward, sideways, vertical);

    guiUpdateCallback(); //We frequently ping the gui update     

    //Flush the OpenCL queue, so GL can use the buffers
    getGLContext()->getCLCommandQueue().finish();
  
    ////////////Lighting shadow map creation////////////////////
    //This stage only needs to be performed once per frame

//    if (_shadowMapping)
//      {
//	_VSMShader.attach();
//
//	//Render each light's shadow map
//	_VSMShader["ProjectionMatrix"] = _light0.getProjectionMatrix();
//	_VSMShader["ViewMatrix"] = _light0.getViewMatrix();	  
//	_light0.shadowFBO().attach();
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	
//	//Enter the render ticks for all objects
//	for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
//	     iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
//	  if ((*iPtr)->shadowCasting() && (*iPtr)->visible())
//	    (*iPtr)->glRender(_light0.shadowFBO(), _light0, RenderObj::SHADOW);
//	
//	_light0.shadowFBO().detach();
//	/////////////MIPMAPPED shadow maps don't seem to work
//	//_light0.shadowTex()->genMipmaps();
//	_light0.shadowTex()->bind(7);
//	
//	_VSMShader.detach();
//      }
    
    ////////3D or Stereo rendering image composition//////////

#ifdef COIL_wiimote
    //Run an update if the wiiMote was connected
    if ((magnet::TrackWiimote::getInstance()).connected())
      {
	Gtk::CheckButton* wiiHeadTrack;
	_refXml->get_widget("wiiHeadTracking", wiiHeadTrack);
	if (wiiHeadTrack->get_active())
	  _camera.setEyeLocation((magnet::TrackWiimote::getInstance()).getHeadPosition());
      }
#endif

    //Bind to the multisample buffer
    if (!_stereoMode)
      {
	drawScene(_renderTarget, _camera);
	_renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight());
      }
    else
      {
	const double eyedist = 6.5;
	Vector eyeDisplacement(0.5 * eyedist, 0, 0);

	Vector currentEyePos = _camera.getEyeLocation();

	Gtk::ComboBox* stereoMode;
	_refXml->get_widget("StereoMode", stereoMode);
	int mode = stereoMode->get_active_row_number();

	switch(mode)
	  {
	  case 0: //Analygraph Red-Cyan
	    _camera.setEyeLocation(currentEyePos + eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    _copyShader["u_Texture0"] = 0;
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	    { 
	      std::tr1::array<GLfloat, 2> arg 
		= {{GLfloat(1) / _camera.getWidth(), 
		    GLfloat(1) / _camera.getHeight()}};
	      _copyShader["u_Scale"] = arg;
	    }
	    glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
	    _copyShader.invoke(); 
	    _copyShader.detach();
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	    _camera.setEyeLocation(currentEyePos - eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.getColorTexture(0)->bind(0);
	    _copyShader.attach();
	    glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE);
	    _copyShader.invoke(); 
	    _copyShader.detach();
	    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	    break;
	  case 1:
	    _camera.setEyeLocation(currentEyePos - eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.blitToScreen(_camera.getWidth() / 2, 
				       _camera.getHeight(), 0, 0, GL_LINEAR);

	    _camera.setEyeLocation(currentEyePos + eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.blitToScreen(_camera.getWidth() / 2, _camera.getHeight(),
				       _camera.getWidth() / 2, 0, GL_LINEAR);	    
	    break;
	  case 2:
	    _camera.setEyeLocation(currentEyePos + eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight()  /2,
				       0, 0, GL_LINEAR);

	    _camera.setEyeLocation(currentEyePos - eyeDisplacement);
	    drawScene(_renderTarget, _camera);
	    _renderTarget.blitToScreen(_camera.getWidth(), _camera.getHeight() / 2,
				       0, _camera.getHeight() / 2, GL_LINEAR);
	    break;
	  default:
	    M_throw() << "Unknown stereo render mode";
	  }
	//Reset the eye position
	_camera.setEyeLocation(currentEyePos);
      }

    getGLContext()->swapBuffers();

    //Check if we're recording and then check that if we're
    //framelocking, check that new data is available
    if ((_snapshot || _record) && (!_simframelock || _newData))
      {
	_newData = false;
	
	std::vector<uint8_t> pixels;
	pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	//Read the pixels into our container
	_renderTarget.getColorTexture()->writeto(pixels);
	
	std::string path;
	{
	  Gtk::FileChooserButton* fileChooser;
	  _refXml->get_widget("snapshotDirectory", fileChooser);
	  path = fileChooser->get_filename();
	}
	
	if (_record || _snapshot)
	  {
	    _snapshot = false;
	    std::ostringstream filename;
	    filename << std::setw(6) <<  std::setfill('0') << std::right << std::dec << _snapshot_counter++;
	    
	    magnet::image::writePNGFile(path + "/" + filename.str() +".png", pixels, 
					_camera.getWidth(), _camera.getHeight(), 4, 1, true, true);
	  }
      }

    ++_frameCounter; 
    _lastFrameTime = _currFrameTime;
    _frameRenderTime = glutGet(GLUT_ELAPSED_TIME) - _currFrameTime;
  }

  void 
  CLGLWindow::drawScene(magnet::GL::FBO& fbo, magnet::GL::Camera& camera)
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    _renderShader.attach();
    _renderShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
    _renderShader["ViewMatrix"] = _camera.getViewMatrix();
    
    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible()) 
	(*iPtr)->glRender(fbo, camera, RenderObj::DEFAULT);

    _renderShader.detach();
    _Gbuffer.detach();
    glDisable(GL_DEPTH_TEST);
    
    ///////////////////////Lighting pass////////////////////////
    //Here we calculate the lighting of every pixel in the scene
    _Gbuffer.getColorTexture(0)->bind(0);
    _Gbuffer.getColorTexture(1)->bind(1);
    _Gbuffer.getColorTexture(2)->bind(2);
    _Gbuffer.getDepthTexture()->bind(3);

    _lightBuffer.attach();
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    //Additive blending of all of the lights contributions
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    _pointLightShader.attach();
    _pointLightShader["colorTex"] = 0;
    _pointLightShader["normalTex"] = 1;
    _pointLightShader["positionTex"] = 2;
    _pointLightShader["depthTex"] = 3;
    _pointLightShader["samples"] = GLint(_samples);
    GLfloat ambient = _ambientIntensity;

    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if (std::tr1::dynamic_pointer_cast<RLight>(*iPtr))
	{
	  std::tr1::shared_ptr<RLight> light 
	    = std::tr1::dynamic_pointer_cast<RLight>(*iPtr);
	  _pointLightShader["ambientLight"] = ambient;
	  _pointLightShader["backColor"] = _backColor;
	  _pointLightShader["lightAttenuation"] = light->getAttenuation();
	  _pointLightShader["lightSpecularExponent"] = light->getSpecularExponent();
	  _pointLightShader["lightSpecularFactor"] = light->getSpecularFactor();
	  _pointLightShader["lightIntensity"] = light->getIntensity();
	  _pointLightShader["lightPosition"] = light->getEyespacePosition(camera);
	  _pointLightShader.invoke();
	  ambient = 0;
	}
    
    _pointLightShader.detach();
    _lightBuffer.detach();
    glDisable(GL_BLEND);

    ///////////////////////Luminance Sampling//////////////////////
    //The light buffer is bound to texture unit 0 for the tone mapping too
    _lightBuffer.getColorTexture()->bind(0);

    _luminanceBuffer.attach();
    _luminanceShader.attach();
    _luminanceShader["colorTex"] = 0;
    _luminanceShader.invoke();
    _luminanceShader.detach();
    _luminanceBuffer.detach();

    //Now we need to generate the mipmaps for the bloom target
    {
      magnet::GL::Texture2D& tex = *_luminanceBuffer.getColorTexture();
      GLsizei currentWidth = tex.getWidth();
      GLsizei currentHeight = tex.getHeight();
      int numLevels = 1 + int(floorf(log2f(fmaxf(currentWidth, currentHeight))));

      //Ensure the luminance buffer is both attached and its color
      //texture bound
      _luminanceBuffer.attach();
      tex.bind(0);

      //Attach the mipmapping shader
      _luminanceMipMapShader.attach();
      _luminanceMipMapShader["luminanceTex"] = 0;
      for (int i=1; i < numLevels; ++i)
	{
	  GLsizei oldWidth = currentWidth;
	  GLsizei oldHeight = currentHeight;
	  //Halve the size of the textures, ensuring they never drop below 1
	  currentWidth /= 2; currentWidth += !currentWidth;
	  currentHeight /= 2; currentHeight += !currentHeight;
	  _glContext->setViewport(0, 0, currentWidth, currentHeight);

	  tex.parameter(GL_TEXTURE_BASE_LEVEL, i - 1);
	  tex.parameter(GL_TEXTURE_MAX_LEVEL, i - 1);
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
				    tex.getGLType(), tex.getGLHandle(), i);

	  //Now generate the mipmap level using a shader
	  std::tr1::array<GLfloat, 2> oldInvDimensions = {{1.0 / oldWidth, 
							   1.0 / oldHeight}};
	  _luminanceMipMapShader["oldInvDimensions"] = oldInvDimensions;
	  std::tr1::array<GLint,2> oldDimensions = {{oldWidth, oldHeight}};
	  _luminanceMipMapShader["oldDimensions"] = oldDimensions;
	  _luminanceMipMapShader.invoke();
	}
      //Rebind mipmap 0 to the framebuffer
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
				tex.getGLType(), tex.getGLHandle(), 0);
      _glContext->setViewport(0, 0, tex.getWidth(), tex.getHeight());

      tex.parameter(GL_TEXTURE_BASE_LEVEL, 0);
      tex.parameter(GL_TEXTURE_MAX_LEVEL, numLevels - 1);
      _luminanceMipMapShader.detach();
      _luminanceBuffer.detach();
    }

    ///////////////////////Tone Mapping///////////////////////////

    fbo.attach();
    _toneMapShader.attach();
    _lightBuffer.getColorTexture()->bind(0);
    _luminanceBuffer.getColorTexture()->bind(1);
    _toneMapShader["color_tex"] = 0;
    _toneMapShader["logLuma"] = 1;
    _toneMapShader["burnout"] = GLfloat((1.0 - _exposure) * _burnoutFactor + _exposure);
    _toneMapShader["exposure"] = GLfloat(_exposure);
    _toneMapShader.invoke();
    _toneMapShader.detach();
    fbo.detach();
    ///////////////////////Forward Shading Pass /////////////////

//    //Enter the forward render ticks for all objects
//    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
//	   = _renderObjsTree._renderObjects.begin();
//	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
//      if ((*iPtr)->visible())
//	(*iPtr)->forwardRender(fbo, camera, , RenderObj::DEFAULT);
//    fbo.detach();
    //////////////////////FILTERING////////////
    //Attempt to perform some filtering

    bool FBOalternate = false;
    magnet::GL::FBO* lastFBO = &fbo;
    
    if (_filterEnable)
      {
       	//Bind the original image to texture unit 0
       	fbo.getColorTexture(0)->bind(0);	
       	//Now bind the texture which has the normals unit 1
       	_Gbuffer.getColorTexture(1)->bind(1);
       	//Positional information is attached to unit 2
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
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    lastFBO->attach();
    _simpleRenderShader.attach();
    _simpleRenderShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
    _simpleRenderShader["ViewMatrix"] = _camera.getViewMatrix();


    //Enter the interface draw for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->interfaceRender(_camera);

    _simpleRenderShader.detach();
    lastFBO->detach();

    glDisable(GL_BLEND);

    //Check if we actually did something and copy the data to the
    //output FBO if needed
    if (lastFBO != &fbo)
      {
	lastFBO->attach();
	fbo.getColorTexture(0)->bind(0);
	glActiveTextureARB(GL_TEXTURE0);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, _camera.getWidth(), 
			    _camera.getHeight());
	lastFBO->detach();
      }

    glEnable(GL_DEPTH_TEST);   
  }

  void CLGLWindow::CallBackReshapeFunc(int w, int h)
  {
    if (!CoilRegister::getCoilInstance().isRunning() || !_readyFlag) return;

    _camera.setHeightWidth(h, w);
    //Update the viewport
    _lightBuffer.resize(w, h);
    _renderTarget.resize(w, h);  
    _Gbuffer.resize(w, h);  
    _filterTarget1.resize(w, h);
    _filterTarget2.resize(w, h);
    _luminanceBuffer.resize(w, h);
    _luminanceBuffer.getColorTexture(0)->genMipmaps();
    std::ostringstream os;
    os << "Coil visualizer (" << w << "," << h << ")";
    setWindowtitle(os.str());
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
	  
	    keyState |= LEFTMOUSE;
	  }
	else
	  keyState &= ~LEFTMOUSE;
	break;
      case GLUT_RIGHT_BUTTON:
	if (state == GLUT_DOWN)
	  {
	    _oldMouseX = x;
	    _oldMouseY = y;
	  
	    keyState |= RIGHTMOUSE;

	    //Now perform a picking selection
	    performPicking(x,y);
	  }
	else
	  keyState &= ~RIGHTMOUSE;
	break;
      case GLUT_MIDDLE_BUTTON:
	if (state == GLUT_DOWN)
	  {
	    _oldMouseX = x;
	    _oldMouseY = y;
	  
	    keyState |= MIDDLEMOUSE;
	  }
	else
	  keyState &= ~MIDDLEMOUSE;
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

    switch (keyState)
      {
      case LEFTMOUSE:
	_camera.mouseMovement(diffX, diffY);
      case RIGHTMOUSE:
      case MIDDLEMOUSE:
      default:
	{}
      }
  
    _oldMouseX = x;
    _oldMouseY = y;
  }

  void 
  CLGLWindow::CallBackKeyboardFunc(unsigned char key, int x, int y)
  {
    keyStates[std::tolower(key)] = true;
  }

  void 
  CLGLWindow::CallBackKeyboardUpFunc(unsigned char key, int x, int y)
  {
    keyStates[std::tolower(key)] = false;
  }

  void
  CLGLWindow::simupdateTick(double t)
  {
    for (;;)
      {
	//Jump out without an update if the window has been killed
	if (!isReady()) return;

	_systemQueue->drainQueue();

	//Block the simulation if _simrun is false or if we're in frame lock
	//and a new frame has not been drawn.
	if (_simrun && (!_simframelock || (_lastUpdateTime != getLastFrameTime()))) break;
      

	//1ms delay to lower CPU usage while blocking, but not to affect framelocked render rates
	timespec sleeptime;
	sleeptime.tv_sec = 0;
	sleeptime.tv_nsec = 1000000;
	nanosleep(&sleeptime, NULL);
      }
    
    //For the updates per second
    ++_updateCounter;

    //Only redraw if the screen has actually refreshed
    if (_lastUpdateTime == getLastFrameTime()) return;
    _lastUpdateTime = getLastFrameTime();

    //Update the simulation data
    {
      magnet::thread::ScopedLock lock(_destroyLock);
      if (!isReady()) return;
      _updateDataSignal();
      _newData = true;

      std::ostringstream os;
      os << "t:" << t;        
      setSimStatus1(os.str());
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
    {
      Gtk::Window* aboutWindow;
      _refXml->get_widget("aboutSplashWindow", aboutWindow);
      aboutWindow->show();
    }

    {
      Gtk::Image* aboutImage;
      _refXml->get_widget("aboutSplashImage", aboutImage);
  
      aboutImage->set(Gdk::Pixbuf::create_from_inline
		      (coilsplash_size, coilsplash));
    }
  }

  void
  CLGLWindow::performPicking(int x, int y)
  {
    _simpleRenderShader.attach();
    _simpleRenderShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
    _simpleRenderShader["ViewMatrix"] = _camera.getViewMatrix();

    _renderTarget.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Flush the OpenCL queue, so GL can use the buffers
    getGLContext()->getCLCommandQueue().finish();
    
    //Perform unique coloring of screen objects
    uint32_t offset = 0;
    //Now render the scene
    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator 
	   iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible())
	(*iPtr)->pickingRender(_renderTarget, _camera, offset);

    unsigned char pixel[4];  
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);  
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel);    
    _renderTarget.detach();

    //Now let the objects know what was picked
    const cl_uint objID = pixel[0] + 256 * (pixel[1] + 256 * (pixel[2] + 256 * pixel[3]));
    offset = 0;
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible())
	(*iPtr)->finishPicking(offset, objID);
  }

  void CLGLWindow::selectRObjCallback() 
  {
    Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
      _renderObjsTree._view->get_selection();

    Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

    Gtk::ScrolledWindow* win;
    _refXml->get_widget("ObjectOptions", win);

    win->remove(); //Clear the current object controls
    if(iter)
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
    switch (_camera.getMode())
      {
      case magnet::GL::Camera::ROTATE_CAMERA:
	_camera.setMode(magnet::GL::Camera::ROTATE_WORLD);
	break;
      case magnet::GL::Camera::ROTATE_WORLD:
	_camera.setMode(magnet::GL::Camera::ROTATE_CAMERA);
	break;
      default:
	M_throw() << "Cannot change camera mode as it's in an unknown mode";
      }
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

      switch (_camera.getMode())
	{
	case magnet::GL::Camera::ROTATE_CAMERA:
	  icon->set(Gdk::Pixbuf::create_from_inline(cammode_fps_size, cammode_fps));
	  break;
	case magnet::GL::Camera::ROTATE_WORLD:
	  icon->set(Gdk::Pixbuf::create_from_inline(cammode_rotate_size, cammode_rotate));
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
      Gtk::Entry* burnoutEntry;
      _refXml->get_widget("BurnoutEntry", burnoutEntry);
      magnet::gtk::forceNumericEntry(*burnoutEntry);
      try {
	_burnoutFactor = boost::lexical_cast<double>(burnoutEntry->get_text());
      } catch(...) {}
    }

    {
      Gtk::Entry* exposureEntry;
      _refXml->get_widget("ExposureEntry", exposureEntry);
      magnet::gtk::forceNumericEntry(*exposureEntry);
      try {
	_exposure = boost::lexical_cast<double>(exposureEntry->get_text());
      } catch(...) {}
    }

    {
      Gtk::Entry* ambientIntensityEntry;
      _refXml->get_widget("AmbientLightIntensity", ambientIntensityEntry);
      magnet::gtk::forceNumericEntry(*ambientIntensityEntry);
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
    
      magnet::gtk::forceNumericEntry(*updateFreq);
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
      if (val.empty()) {val = "50"; simunits->set_text("50"); }
      _camera.setSimUnitLength(boost::lexical_cast<double>(val));
    }

    {
      Gtk::Entry* pixelPitch;
      _refXml->get_widget("pixelPitch", pixelPitch);
      std::string val = pixelPitch->get_text();
      if (val.empty()) {val = "0.25"; pixelPitch->set_text("0.25"); }
      _camera.setPixelPitch(boost::lexical_cast<double>(val) / 10);
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
  CLGLWindow::HeadReset()
  {
    _camera.setEyeLocation(Vector(0, 0, _camera.getEyeLocation()[2]));
    _camera.setFOVY(60.f, false);
  }
}

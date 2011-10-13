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

extern const guint8 coilsplash[];
extern const size_t coilsplash_size;

//The glade xml file is "linked" into a binary file and stuffed in the executable, these are the symbols to its data
extern const char _binary_clwingtk_gladexml_start[];
extern const char _binary_clwingtk_gladexml_end[];
extern const guint8 coilicon[];
extern const size_t coilicon_size;


namespace coil {
  CLGLWindow::CLGLWindow(std::string title,
			 double updateIntervalValue,
			 bool dynamo
			 ):
    _systemQueue(new magnet::thread::TaskQueue),
    _updateIntervalValue(updateIntervalValue),
    _glContext(NULL),
    keyState(DEFAULT),
    windowTitle(title),
    _frameCounter(0),
    _updateCounter(0),
    _mouseSensitivity(0.3),
    _moveSensitivity(0.001),
    _specialKeys(0),
    _shadowMapping(false),
    _shadowIntensity(0.8),
    _simrun(false),
    _simframelock(false),
    _snapshot(false),
    _record(false),
    _PNGFileFormat(true),
    _fpsLimit(true),
    _fpsLimitValue(35),
    _filterEnable(true),
    _analygraphMode(false),
    _snapshot_counter(0),
    _dynamo(dynamo)
  {
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

    {//////////////Glade XML loader 
      Glib::ustring glade_data
	(reinterpret_cast<const char *>(_binary_clwingtk_gladexml_start), 
	 _binary_clwingtk_gladexml_end
	 -_binary_clwingtk_gladexml_start);
    
      _refXml = Gtk::Builder::create_from_string(glade_data);
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

    {//////Place light button
      Gtk::Button* lightButton;    
      _refXml->get_widget("lightLocation", lightButton); 

      lightButton->signal_clicked()
	.connect(sigc::mem_fun(*this, &CLGLWindow::lightPlaceCallback));
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
      Gtk::SpinButton* updateButton;
      _refXml->get_widget("updateFreq", updateButton);
      updateButton->set_value(_updateIntervalValue);
      updateButton->signal_value_changed()
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


    {///////Light FOV setting
      Gtk::HScale* FOVscale;
      _refXml->get_widget("lightFOVScale", FOVscale);
      FOVscale->set_value(_light0.getFOVY());
      FOVscale->signal_value_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
    }

    ///////////////////////Render Pipeline//////////////////////////////////
    {
      ///////////////////////Shadow Mapping//////////////////////////////////
      {
	Gtk::CheckButton* shadowmapEnable;
	_refXml->get_widget("shadowmapEnable", shadowmapEnable);

	shadowmapEnable->set_active(_shadowMapping);
	shadowmapEnable->signal_toggled()
	  .connect(sigc::mem_fun(this, &CLGLWindow::shadowEnableCallback));
      }
    
      {
	Gtk::SpinButton* shadowmapSize;
	_refXml->get_widget("shadowmapSize", shadowmapSize);
	shadowmapSize->set_value(1024);
	shadowmapSize->signal_value_changed()
	  .connect(sigc::mem_fun(this, &CLGLWindow::shadowEnableCallback));
      }
    
      {//Setup the shadow intensity
	Gtk::VolumeButton* shadowButton;
	_refXml->get_widget("shadowIntensity", shadowButton);
	shadowButton->set_value(_shadowIntensity);
      
	shadowButton->signal_value_changed()
	  .connect(sigc::mem_fun(this, &CLGLWindow::shadowIntensityCallback));
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
	  Gtk::CheckButton* analygraphEnable;
	  _refXml->get_widget("analygraphMode", analygraphEnable);
	  analygraphEnable->signal_toggled()
	    .connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
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
  CLGLWindow::shadowEnableCallback()
  {
    Gtk::CheckButton* shadowmapEnable;
    _refXml->get_widget("shadowmapEnable", shadowmapEnable);
  
    _shadowMapping = shadowmapEnable->get_active();


    if (_shadowMapping)
      {
	Gtk::SpinButton* shadowmapSize;
	_refXml->get_widget("shadowmapSize", shadowmapSize);
      
	_light0.shadowFBO().resize(shadowmapSize->get_value(), shadowmapSize->get_value());
      }
  }

  void
  CLGLWindow::init()
  {
    magnet::thread::ScopedLock lock(_destroyLock);

    if (_readyFlag) return;

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

    //  //Test volume render object
    
//    std::tr1::shared_ptr<RVolume> vol(new RVolume("Test Volume"));
//    _renderObjsTree._renderObjects.push_back(vol);

//    glutInitContextVersion(3, 3);
//    
//#ifdef MAGNET_DEBUG
//    glutInitContextFlags(GLUT_CORE_PROFILE | GLUT_DEBUG);
//#else
//    glutInitContextFlags(GLUT_CORE_PROFILE);
//#endif
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(0, 0);

    CoilRegister::getCoilInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

    _glContext = &magnet::GL::Context::getContext();

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    //Blend colors using the alpha channel
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

    //Switch on line aliasing
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    //Setup the viewport
    CallBackReshapeFunc(800, 600);

    _light0 = magnet::GL::Light(Vector(0.8f,  1.5f, 0.8f),//Position
				Vector(0.0f, 0.0f, 0.0f),//Lookat
				75.0f//Beam angle
				);
  
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

    _light0.init();
    
    _renderShader.build();
    _deferredShader.build();
    _VSMShader.build();
    _simpleRenderShader.build();

    { 
      //Build render buffer
      std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
      depthTexture->init(_camera.getWidth(), _camera.getHeight(), GL_DEPTH_COMPONENT);
      depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
      
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGB);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      _renderTarget.init();
      _renderTarget.attachTexture(colorTexture, 0);
      _renderTarget.attachTexture(depthTexture);
    }

    {
      //Build G buffer
      std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
      depthTexture->init(_camera.getWidth(), _camera.getHeight(), GL_DEPTH_COMPONENT);
      depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
      
      std::tr1::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D);
      colorTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
      colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      
      std::tr1::shared_ptr<magnet::GL::Texture2D> normalTexture(new magnet::GL::Texture2D);
      normalTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
      normalTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      normalTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      normalTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      normalTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

      std::tr1::shared_ptr<magnet::GL::Texture2D> posTexture(new magnet::GL::Texture2D);
      posTexture->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
      posTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      posTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      posTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      posTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
      _Gbuffer.init();
      _Gbuffer.attachTexture(colorTexture, 0);
      _Gbuffer.attachTexture(normalTexture, 1);
      _Gbuffer.attachTexture(posTexture, 2);
      _Gbuffer.attachTexture(depthTexture);
    }

    //Now init the render objects  
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->init(_systemQueue);
  
    Console& _console = static_cast<Console&>(*_renderObjsTree._renderObjects[_consoleID]);
    _console << "Welcome to coil, part of the dynamo simulator..." 
	     << Console::end();

    initGTK();

    //  //Fabian Test
    //  vol->loadRawFile("/home/mjki2mb2/Desktop/Output.raw", 300, 300, 300, 1);
    //  
    //bonsai plant test
    //  vol->loadRawFile("bonsai.raw", 256, 256, 256, 1);
    //
    //  //Cadaver
    //  vol->loadRawFile("cadaver512x512x106.raw", 512, 512, 106, 2);
    //
    //  //Male 
    //  vol->loadRawFile("Male128x256x256.raw", 128, 256, 256, 1);
    //
    //  //Female
    //  vol->loadRawFile("female384x384x240.raw", 384, 384, 240, 1);

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

    /////////////////OpenCL

    getGLContext().getCLCommandQueue().finish();

    ///////////////////OpenGL
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->deinit();

    _renderObjsTree._renderObjects.clear();

    _renderTarget.deinit();
    _Gbuffer.deinit();
    _filterTarget1.deinit();
    _filterTarget2.deinit();
    _renderShader.deinit();
    _deferredShader.deinit();	
    _VSMShader.deinit();
    _simpleRenderShader.build();

    _light0.deinit();
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
    glFinish();//Finish with the GL buffers

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

#ifdef COIL_wiimote
    //Run an update if the wiiMote was connected
    if ((magnet::TrackWiimote::getInstance()).connected())
      {
	Gtk::CheckButton* wiiHeadTrack;
	_refXml->get_widget("wiiHeadTracking", wiiHeadTrack);
	if (wiiHeadTrack->get_active())
	  _camera.setHeadLocation((magnet::TrackWiimote::getInstance()).getHeadPosition());
      }
#endif

    //Flush the OpenCL queue, so GL can use the buffers
    getGLContext().getCLCommandQueue().finish();
  
    //Prepare for the GL render
    if (_shadowMapping)
      {
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_BLEND);
	_VSMShader.attach();

	//Render each light's shadow map
	_VSMShader["ProjectionMatrix"] = _light0.getProjectionMatrix();
	_VSMShader["ViewMatrix"] = _light0.getViewMatrix();	  
	_light0.shadowFBO().attach();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	
	//Enter the render ticks for all objects
	for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	     iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
	  if ((*iPtr)->shadowCasting() && (*iPtr)->visible())
	    (*iPtr)->glRender(_light0.shadowFBO(), _light0, RenderObj::SHADOW);
	
	_light0.shadowFBO().detach();
	_light0.shadowTex().genMipmaps();
	_light0.shadowTex().bind(7);
	
	_VSMShader.detach();
	glEnable(GL_BLEND);
	glEnable(GL_ALPHA_TEST);
      }
      
    _renderTarget.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    _renderTarget.detach();
    
    //Bind to the multisample buffer
    if (_analygraphMode)
      {
	const double eyedist = 6.5;
	Vector eyeDisplacement(0.5 * eyedist, 0, 0);
	  
	glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
	drawScene(_renderTarget, _camera, eyeDisplacement);

	glClear(GL_DEPTH_BUFFER_BIT);
	glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE);
	drawScene(_renderTarget, _camera, -eyeDisplacement);

	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
      }
    else
      drawScene(_renderTarget, _camera, Vector(0,0,0));

    //////////////FILTERING////////////
    bool FBOalternate = false;

    magnet::GL::FBO* lastFBO = &_renderTarget;
    if (_filterEnable && !_filterStore->children().empty())
      {
	glDisable(GL_DEPTH_TEST);

	//Bind the original image to texture (unit 0)
	_renderTarget.getColorTexture(0).bind(0);	
	//Now bind the texture which has the normals (unit 1)
	_Gbuffer.getColorTexture(1).bind(1);
	//High quality depth information is attached to (unit 2)
	_Gbuffer.getDepthTexture().bind(2);

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
		lastFBO->getColorTexture().bind(3);
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
	glEnable(GL_DEPTH_TEST);
      }

    //Now blit the stored scene to the screen
    lastFBO->blitToScreen(_camera.getWidth(), _camera.getHeight());

    //We clear the depth as merely disabling gives artifacts
    glClear(GL_DEPTH_BUFFER_BIT); 

    _simpleRenderShader.attach();
    _simpleRenderShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
    _simpleRenderShader["ViewMatrix"] = _camera.getViewMatrix();

    //Enter the interface draw for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      (*iPtr)->interfaceRender(_camera);

    _simpleRenderShader.detach();

    getGLContext().swapBuffers();

    //Check if we're recording and then check that if we're
    //framelocking, check that new data is available
    if ((_snapshot || _record) && (!_simframelock || _newData))
      {
	_newData = false;

	std::vector<uint8_t> pixels;
	pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	//Read the pixels into our container
	lastFBO->getColorTexture().getTexture(pixels);

	std::string path;
	{
	  Gtk::FileChooserButton* fileChooser;
	  _refXml->get_widget("snapshotDirectory", fileChooser);
	  path = fileChooser->get_filename();
	}

	if (_snapshot)
	  {
	    _snapshot = false;
	    magnet::image::writePNGFile(path + "/snapshot.png", pixels, _camera.getWidth(), 
					_camera.getHeight(), 4, 9, false, true);
	  }
	
	if (_record)
	  {
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
  CLGLWindow::drawScene(magnet::GL::FBO& fbo, magnet::GL::Camera& camera, Vector eyeDisplacement)
  {
    _Gbuffer.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    _renderShader.attach();
    _renderShader["ProjectionMatrix"] = _camera.getProjectionMatrix(eyeDisplacement);
    _renderShader["ViewMatrix"] = _camera.getViewMatrix(eyeDisplacement);

    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible()) (*iPtr)->glRender(fbo, camera, RenderObj::DEFAULT);

    _renderShader.detach();
    _Gbuffer.detach();

    _renderTarget.attach();
    
    _deferredShader.attach();
    _Gbuffer.getColorTexture(0).bind(0);
    _deferredShader["colorTex"] = 0;
    _Gbuffer.getColorTexture(1).bind(1);
    _deferredShader["normalTex"] = 1;
    _Gbuffer.getColorTexture(2).bind(2);
    _deferredShader["positionTex"] = 2;

    {
      magnet::math::Vector vec = _light0.getEyeLocation();
      std::tr1::array<GLfloat, 4> lightPos = {{vec[0], vec[1], vec[2], 1.0}};
      std::tr1::array<GLfloat, 4> lightPos_eyespace 
	= camera.getViewMatrix(eyeDisplacement) * lightPos;
      magnet::math::Vector vec2(lightPos_eyespace[0], lightPos_eyespace[1], lightPos_eyespace[2]);
      _deferredShader["lightPosition"] = vec2;
    }

    _deferredShader["ShadowMap"] = 7;
    _deferredShader["ShadowIntensity"] = _shadowIntensity;
    _deferredShader["ShadowMapping"] = _shadowMapping;
    if (_shadowMapping)
      _deferredShader["ShadowMatrix"] = _light0.getShadowTextureMatrix() * camera.getViewMatrix(eyeDisplacement).inverse();

    _deferredShader.invoke();

    _deferredShader.detach();
    _renderTarget.detach();
  }

  void CLGLWindow::CallBackReshapeFunc(int w, int h)
  {
    if (!CoilRegister::getCoilInstance().isRunning() || !_readyFlag) return;

    _camera.setHeightWidth(h, w);
    //Update the viewport
    _renderTarget.resize(w, h);  
    _Gbuffer.resize(w, h);  
    _filterTarget1.resize(w, h);
    _filterTarget2.resize(w, h);
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
    ++_updateCounter;//For the updates per second

    for (;;)
      {
	_systemQueue->drainQueue();

	//Block the simulation if _simrun is false or if we're in frame lock
	//and a new frame has not been drawn.
	if (_simrun && (!_simframelock || (_lastUpdateTime != getLastFrameTime()))) break;
      
	//Jump out without an update if the window has been killed
	if (!isReady()) return;

	//1ms delay to lower CPU usage while blocking, but not to affect framelocked render rates
	timespec sleeptime;
	sleeptime.tv_sec = 0;
	sleeptime.tv_nsec = 1000000;
	nanosleep(&sleeptime, NULL);
      }

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
  CLGLWindow::lightPlaceCallback()
  { _light0 = _camera; }

  void 
  CLGLWindow::shadowIntensityCallback(double val)
  {
    _shadowIntensity = val;
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
    _renderTimeout = Glib::signal_timeout().connect(sigc::mem_fun(this, &CLGLWindow::CallBackIdleFunc), 
						    _fpsLimit ? 1000 / _fpsLimitValue : 10, 
						    Glib::PRIORITY_DEFAULT_IDLE);
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
    //We need a non-multisampled FBO, just use one of the filter FBO's
    _filterTarget1.attach();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_DITHER);
    glShadeModel(GL_FLAT);
    
    //Flush the OpenCL queue, so GL can use the buffers
    getGLContext().getCLCommandQueue().finish();
    
    //Perform unique coloring of screen objects
    uint32_t offset = 0;
    //Now render the scene
    //Enter the render ticks for all objects
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator 
	   iPtr = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible())
	(*iPtr)->pickingRender(_filterTarget1, _camera, offset);

    unsigned char pixel[4];  
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);  
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel);
    
    std::cout << "Pixel read gave " 
	      << int(pixel[0]) << " "
	      << int(pixel[1]) << " "
	      << int(pixel[2]) << " "
	      << int(pixel[3]) << std::endl;

    _filterTarget1.detach();
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glEnable(GL_DITHER);
    glShadeModel(GL_SMOOTH);

    //Now let the objects know what was picked
    const cl_uint objID = pixel[0] + 256 * (pixel[1] + 256 * (pixel[2] + 256 * pixel[3]));
    offset = 0;
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr 
	   = _renderObjsTree._renderObjects.begin();
	 iPtr != _renderObjsTree._renderObjects.end(); ++iPtr)
      if ((*iPtr)->visible())
	(*iPtr)->finishPicking(offset, objID);

    //Debugging output for the picking buffer
    //_filterTarget1.blitToScreen(_camera.getWidth(), _camera.getHeight());
    //getGLContext().swapBuffers();
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
  CLGLWindow::guiUpdateCallback()
  {
    {///////light FOV setting
      Gtk::HScale* FOVscale;
      _refXml->get_widget("lightFOVScale", FOVscale);
      _light0.setFOVY(FOVscale->get_value());
    }

    {//Dynamo particle sync checkbox
      Gtk::CheckButton* btn;
      _refXml->get_widget("forceParticleSync", btn);
    
      _particleSync = btn->get_active();
    }

    {//Filter enable/disable
      Gtk::CheckButton* btn;
      _refXml->get_widget("filterEnable", btn);
    
      _filterEnable = btn->get_active();
    }

    {//Sim Update Frequency Control
      Gtk::SpinButton* updateButton;
      _refXml->get_widget("updateFreq", updateButton);
    
      if (updateButton->get_value() <= 0)
	updateButton->set_value(0.000001);
    
      _updateIntervalValue = updateButton->get_value();
    }

    {//Analygraph work
      Gtk::CheckButton* btn;
      _refXml->get_widget("analygraphMode", btn);    
      _analygraphMode = btn->get_active();
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
      os << _camera.getHeadLocation()[0] << "cm";
      XHead->set_text(os.str());
      os.str("");
      os << _camera.getHeadLocation()[1] << "cm";
      YHead->set_text(os.str());
      os.str("");
      os << _camera.getHeadLocation()[2] << "cm";
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
  CLGLWindow::HeadReset()
  {
    _camera.setHeadLocation(Vector(0,0,_camera.getHeadLocation()[2]));
    _camera.setFOVY(60.f, false);
  }
}

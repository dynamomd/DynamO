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
#include <coil/RenderObj/Cylinders.hpp>

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
  _shaderPipeline(false),
  _shadowMapping(true),
  _shadowIntensity(0.8),
  _simrun(false),
  _simframelock(false),
  _snapshot(false),
  _record(false),
  _showLight(false),
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

CLGLWindow::~CLGLWindow()
{}

void 
CLGLWindow::initOpenGL()
{
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowSize(800, 600);
  glutInitWindowPosition(0, 0);

  CoilRegister::getCoilInstance().CallGlutCreateWindow(windowTitle.c_str(), this);

  _glContext = &magnet::GL::Context::getContext();

  glViewport(0, 0, 800, 600);

  if (!GLEW_VERSION_2_0)
    std::runtime_error("OpenGL 2.0 is not supported, coil is not supported in this mode");

  if (!GLEW_EXT_framebuffer_object)
    std::runtime_error("Frame buffers are not supported, coil is not supported without them");

  if (!GLEW_ARB_vertex_buffer_object)
    std::runtime_error("Vertex Buffer Objects are not supported by your GPU/Driver, this is critical for Coil."); 

  //Check for all of the shader pipeline support
  _shaderPipeline = true;
  if (!GLEW_ARB_depth_texture || !GLEW_ARB_shadow)
    {
      std::cout << "GL_ARB_depth_texture or GL_ARB_shadow not supported.\n";
      _shaderPipeline = false;
    }

  else if (!GLEW_ARB_fragment_program || !GLEW_ARB_vertex_program
	   || !GLEW_ARB_fragment_shader || !GLEW_ARB_vertex_shader)
    {
      std::cout << "OpenGL driver doesn't support programmable shaders.\n";
      _shaderPipeline = false;
    }
  
  if (!_shaderPipeline)
    std::cout << "Shader pipeline disabled.\n"
	      << "This also disables all nice effects and alot of speed.\n";

  
  _pickingEnabled = true;
  {
    GLint bits;
    glGetIntegerv(GL_ALPHA_BITS, &bits);
    if (bits < 8) _pickingEnabled = false;
    glGetIntegerv(GL_RED_BITS, &bits);
    if (bits < 8) _pickingEnabled = false;
    glGetIntegerv(GL_GREEN_BITS, &bits);
    if (bits < 8) _pickingEnabled = false;
    glGetIntegerv(GL_BLUE_BITS, &bits);
    if (bits < 8) _pickingEnabled = false;
  }

  if (!_pickingEnabled)
    std::cout << "\nPicking won't work! Your screen color depth is too low! We need/want 32bits (24bit color plus 8bit alpha)";

  _viewPortInfo->buildMatrices();

  glDrawBuffer(GL_BACK);

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);

  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  glEnable(GL_BLEND);
  //Blend colors using the alpha channel
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

  //Make OpenGL renormalize lighting vectors for us (incase we use glScale)
  glEnable(GL_NORMALIZE);

  //Switch on line aliasing
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  //glEnable(GL_CULL_FACE);
  //glCullFace(GL_BACK);
  //glFrontFace(GL_CCW); //The default

  //Both the front and back materials track the current color
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL); //and enable it
  
  glShadeModel(GL_SMOOTH);

  //Setup the viewport
  CallBackReshapeFunc(800, 600);

  glReadBuffer(GL_BACK);
  glPixelStorei(GL_PACK_ALIGNMENT, 4);

  //Light our scene!
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  //Light number one
  //Position is set in the CameraSetup!
  //Ambient lighting
  GLfloat ambient_light[] = {0.0f, 0.0f, 0.0f, 1.0f}; 
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);

  //Light0 parameters
  _light0.reset(new magnet::GL::LightInfo(Vector(0.8f,  1.5f, 0.8f),//Position
					  Vector(0.0f, 0.0f, 0.0f),//Lookat
					  GL_LIGHT0,
					  75.0f//Beam angle
					  ));
  
  _light0->buildMatrices();

  glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0f);
  glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0f);
  glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.1f);
  GLfloat ambientLight0[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight0);

  //Default material parameters
  GLfloat specReflection[] = { 1.0f, 1.0f, 1.0f, 1.0f };
  glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
  glMateriali(GL_FRONT, GL_SHININESS, 25);


  //GLfloat specular[] = {1.0, 0.0, 0.0, 1.0};
  //glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

  //Setup the keyboard controls
  glutIgnoreKeyRepeat(1);

  _lastUpdateTime = _lastFrameTime = _FPStime = glutGet(GLUT_ELAPSED_TIME);
  _frameRenderTime = 0;
  //Build the offscreen rendering FBO's
  _renderTarget.reset(new magnet::GL::FBO());
  _renderTarget->init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());

  if (_shaderPipeline)
    {
      _filterTarget1.init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
      _filterTarget2.init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
      _normalsFBO.init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight(), GL_RGBA);
      _shadowFBO.init(1024);
      _renderShader.build();
      _depthRenderShader.build();
      _nrmlShader.build();
    }

  //Now init the render objects  
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initOpenGL();
  
  coil::Console& _console 
    = static_cast<coil::Console&>(*RenderObjects[_consoleID]);
  _console << "Welcome to coil, part of the dynamo simulator..." 
	   << coil::Console::end();
}

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
CLGLWindow::initOpenCL()
{
  //Force the OpenCL platform to be constructed
  getGLContext().getCLPlatform();
  
  //Now init the render objects  
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initOpenCL();
}

//The glade xml file is "linked" into a binary file and stuffed in the executable, these are the symbols to its data
extern const char _binary_clwingtk_gladexml_start[];
extern const char _binary_clwingtk_gladexml_end[];
extern const guint8 coilicon[];
extern const size_t coilicon_size;

void
CLGLWindow::initGTK()
{
  _filterModelColumns.reset(new FilterModelColumnsType);
  _renderObjModelColumns.reset(new RenderObjModelColumnsType);

  {////////Glade XML loader 
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

  {//////Show light checkbox
    Gtk::CheckButton* lightShowButton;    
    _refXml->get_widget("lightShow", lightShowButton); 

    lightShowButton->signal_toggled()
      .connect(sigc::mem_fun(*this, &CLGLWindow::lightShowCallback));
  }

  {//////Render the surface normals checkbox
    Gtk::CheckButton* Button;    
    _refXml->get_widget("showNormals", Button); 

    Button->signal_toggled()
      .connect(sigc::mem_fun(*this, &CLGLWindow::renderNormalsCallback));
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

  {///////File format selection
    Gtk::RadioButton* radioButton;
    _refXml->get_widget("snapshotBMP", radioButton);
    radioButton->set_active(false);
    radioButton->signal_toggled()
      .connect(sigc::mem_fun(this, &CLGLWindow::snapshotFileFormatCallback));
    _refXml->get_widget("snapshotPNG", radioButton);
    radioButton->set_active(true);
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
    FOVscale->set_value(_light0->getFOVY());
    FOVscale->signal_value_changed()
      .connect(sigc::mem_fun(this, &CLGLWindow::guiUpdateCallback));
  }

  ///////////////////////Render Pipeline//////////////////////////////////
  if (_shaderPipeline)
    {

      {//Enable the whole shader frame
	Gtk::Frame* shaderFrame;
	_refXml->get_widget("RenderPipelineFrame", shaderFrame);
	shaderFrame->set_sensitive(true);
      }
      
      {//Setup the shader enable
	Gtk::CheckButton* shaderEnable;
	_refXml->get_widget("ShaderPipelineEnable", shaderEnable);
	
	shaderEnable->set_active(true);

	shaderEnable->signal_toggled().connect(sigc::mem_fun(this, &CLGLWindow::pipelineEnableCallback));
      }

      ///////////////////////Multisampling (anti-aliasing)//////////////////////////////////
      GLint maxSamples = magnet::GL::MultisampledFBO::getSupportedSamples();

      if (maxSamples > 1)
	{//Offer anti aliasing
	  {//Turn on the antialiasing box
	    Gtk::HBox* multisampleBox;
	    _refXml->get_widget("multisampleBox", multisampleBox);
	    multisampleBox->set_sensitive(true);
	  }

	  //Connect the anti aliasing checkbox
	  Gtk::CheckButton* multisampleEnable;
	  _refXml->get_widget("multisampleEnable", multisampleEnable);
	  multisampleEnable->signal_toggled()
	    .connect(sigc::mem_fun(*this, &CLGLWindow::multisampleEnableCallback));
	  
	  
	  Gtk::ComboBox* aliasSelections;
	  _refXml->get_widget("multisampleLevels", aliasSelections);

	  struct aliasColumns : public Gtk::TreeModel::ColumnRecord
	  {
	    aliasColumns() { add(m_col_id); }
	    Gtk::TreeModelColumn<int> m_col_id;
	  };

	  aliasColumns vals;
	  Glib::RefPtr<Gtk::ListStore> m_refTreeModel = Gtk::ListStore::create(vals);
	  aliasSelections->set_model(m_refTreeModel);

	  Gtk::TreeModel::Row row;
	  int lastrow = -1;
	  GLint currentSamples = maxSamples;
	  for ( ; currentSamples > 1; currentSamples >>= 1)
	    {
	      row = *(m_refTreeModel->prepend());
	      row[vals.m_col_id] = currentSamples;
	      ++lastrow;
	    }

	  aliasSelections->pack_start(vals.m_col_id);

	  //Activate a multisample of 2<<(2)=8 by default
	  aliasSelections->set_active(std::min(lastrow, 2));

	  multisampleEnable->set_active(true);
	  
	  _renderTarget.reset(new magnet::GL::MultisampledFBO
			      (2 << aliasSelections->get_active_row_number()));
	  
	  _renderTarget->init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());

	  aliasSelections->signal_changed()
	    .connect(sigc::mem_fun(*this, &CLGLWindow::multisampleEnableCallback));
	}

      ///////////////////////Shadow Mapping//////////////////////////////////
      {
	Gtk::CheckButton* shadowmapEnable;
	_refXml->get_widget("shadowmapEnable", shadowmapEnable);
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
	  coil::filter::populateComboBox(filterSelectBox);
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
	  os << _viewPortInfo->getSimUnitLength();
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
	  os << _viewPortInfo->getPixelPitch() * 10;
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
    
    //Build the store
    _renderObjStore = Gtk::ListStore::create(*_renderObjModelColumns);
    
    //Setup the filter store
    _refXml->get_widget("renderObjView", _renderObjView);
    _renderObjView->set_model(_renderObjStore);
    _renderObjView->append_column("Visible", _renderObjModelColumns->m_visible);
    _renderObjView->append_column("Object Name", _renderObjModelColumns->m_name);
    
    //Populate the render object view
    rebuildRenderView();
    selectRObjCallback();

    //////Connect the view to the select callback
    {
      Glib::RefPtr<Gtk::TreeSelection> treeSelection
	= _renderObjView->get_selection();
      
      treeSelection->signal_changed()
	.connect(sigc::mem_fun(this, &CLGLWindow::selectRObjCallback));
    }
    
    {///Connect the control buttons
      Gtk::Button* btn;
      _refXml->get_widget("robjDelete", btn);
      btn->signal_clicked()
	.connect(sigc::mem_fun(this, &CLGLWindow::deleteRObjCallback));

      _refXml->get_widget("robjEdit", btn);
      btn->signal_clicked()
	.connect(sigc::mem_fun(this, &CLGLWindow::editRObjCallback));

      _refXml->get_widget("robjAdd", btn);
      btn->signal_clicked()
	.connect(sigc::mem_fun(this, &CLGLWindow::addRObjCallback));
    }
    {
      Gtk::ToggleButton* btn;
      _refXml->get_widget("robjVisible", btn);
      btn->signal_toggled()
	.connect(sigc::mem_fun(this, &CLGLWindow::visibleRObjCallback));
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

  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initGTK();
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
CLGLWindow::pipelineEnableCallback()
{
  {
    Gtk::CheckButton* shaderEnable;
    _refXml->get_widget("ShaderPipelineEnable", shaderEnable);  
    _shaderPipeline = shaderEnable->get_active();
  }

  {
    Gtk::VBox* shaderPipelineOptions;
    _refXml->get_widget("shaderPipelineOptions", shaderPipelineOptions);
    shaderPipelineOptions ->set_sensitive(_shaderPipeline); 
  }
}

void 
CLGLWindow::multisampleEnableCallback()
{
  Gtk::CheckButton* multisampleEnable;
  _refXml->get_widget("multisampleEnable", multisampleEnable);
  if (multisampleEnable->get_active())
    {
      Gtk::ComboBox* aliasSelections;
      _refXml->get_widget("multisampleLevels", aliasSelections);

      _renderTarget.reset(new magnet::GL::MultisampledFBO(2 << aliasSelections->get_active_row_number()));
      _renderTarget->init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
    }
  else
    {
      _renderTarget.reset(new magnet::GL::FBO());
      _renderTarget->init(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
    }
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
      
      _shadowFBO.resize(shadowmapSize->get_value(), shadowmapSize->get_value());
    }
}

void
CLGLWindow::init()
{
  magnet::thread::ScopedLock lock(_destroyLock);

  if (_readyFlag) return;

  //First render object is the ground
  RenderObjects.push_back(new RFunction((size_t)64,
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


  //Second render object is the console
  _consoleID = RenderObjects.size();
  std::tr1::array<GLfloat, 3> textcolor  = {{0.5, 0.5, 0.5}};
  RenderObjects.push_back(new coil::Console(textcolor));

//  //Test volume render object
//  RenderObjects.push_back(new coil::RVolume("Test Volume"));

  //Test cylinder render object
  RenderObjects.push_back(new coil::RCylinders(100000, "Cylinder Test Objects"));

  _viewPortInfo 
    = magnet::thread::RefPtr<magnet::GL::ViewPort>(new magnet::GL::ViewPort);

  //Inform objects about the accessory objects, like the console or the pointer
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    //This is the console and the simulation/system task queue
    (*iPtr)->accessoryData(RenderObjects[_consoleID], _systemQueue, _viewPortInfo); 
  
  initOpenGL();
  initOpenCL();
  initGTK();

//  //Fabian Test
//  RenderObjects.back().as<coil::RVolume>()
//    .loadRawFile("/home/mjki2mb2/Desktop/Output.raw", 300, 300, 300, 1);
//  
//  //bonsai plant test
//  RenderObjects.back().as<coil::RVolume>()
//    .loadRawFile("bonsai.raw", 256, 256, 256, 1);
//
//  //Cadaver
//  RenderObjects.back().as<coil::RVolume>()
//    .loadRawFile("cadaver512x512x106.raw", 512, 512, 106, 2);
//
//  //Male 
//  RenderObjects.back().as<coil::RVolume>()
//    .loadRawFile("Male128x256x256.raw", 128, 256, 256, 1);
//
//  //Female
//  RenderObjects.back().as<coil::RVolume>()
//    .loadRawFile("female384x384x240.raw", 384, 384, 240, 1);

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
      delete static_cast<coil::filter*>(tmp_ptr);
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
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->releaseCLGLResources();

  RenderObjects.clear();

  _renderTarget->deinit();
  _filterTarget1.deinit();
  _filterTarget2.deinit();
  _normalsFBO.deinit();
  _shadowFBO.deinit();
  _renderShader.deinit();
  _depthRenderShader.deinit();
  _nrmlShader.deinit();

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

//  const float speed = 1000;
//  _light0 = lightInfo(Vector(1.5f*std::cos(_currFrameTime/speed), 1.5f, 
//			     1.5f * std::sin(_currFrameTime/speed)), 
//		      Vector(0.0f, 0.0f, 0.0f), GL_LIGHT0);

  //Run every objects OpenCL stage
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->clTick();

  //Camera Positioning

  float moveAmp  = (_currFrameTime - _lastFrameTime) * _moveSensitivity;      
  float forward  = moveAmp * ( keyStates[static_cast<size_t>('w')] 
			       - keyStates[static_cast<size_t>('s')]);
  float sideways = moveAmp * ( keyStates[static_cast<size_t>('d')] 
			       - keyStates[static_cast<size_t>('a')]);
  float vertical = moveAmp * ( keyStates[static_cast<size_t>('q')] 
			       - keyStates[static_cast<size_t>('z')]);
  _viewPortInfo->CameraUpdate(forward, sideways, vertical);

  guiUpdateCallback(); //We frequently ping the gui update     

#ifdef COIL_wiimote
  //Run an update if the wiiMote was connected
  if ((magnet::TrackWiimote::getInstance()).connected())
    {
      {
	Gtk::CheckButton* wiiHeadTrack;
	_refXml->get_widget("wiiHeadTracking", wiiHeadTrack);
	if (wiiHeadTrack->get_active())
	  _viewPortInfo->setHeadLocation((magnet::TrackWiimote::getInstance()).getHeadPosition());
      }
    }
#endif

  //Flush the OpenCL queue, so GL can use the buffers
  getGLContext().getCLCommandQueue().finish();
  
  //Prepare for the GL render
  if (_shaderPipeline)
    {
      if (_shadowMapping)
	{
	  _depthRenderShader.attach();
	  //////////////////Pass 1//////////////////
	  ///Here we draw from the lights perspective
	  _light0->loadMatrices();
	  
	  //Setup the FBO for shadow maps
	  _shadowFBO.setup();
	  
#ifdef GL_VERSION_1_1
	  glEnable(GL_POLYGON_OFFSET_FILL);
	  glPolygonOffset(5.0f, 10.0f);
#endif 

	  drawScene(_shadowFBO);

#ifdef GL_VERSION_1_1
	  glDisable(GL_POLYGON_OFFSET_FILL);
#endif
	  _shadowFBO.restore();
	  _shadowFBO.getDepthTexture().bind(7);
	}
      
      //Bind to the multisample buffer
      _renderTarget->attach();
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

      _renderShader["ShadowMap"] = 7;
      _renderShader["ShadowIntensity"] = _shadowIntensity;
      _renderShader["xPixelOffset"] = 1.0f / _viewPortInfo->getWidth();
      _renderShader["yPixelOffset"] = 1.0f / _viewPortInfo->getHeight();
      _renderShader["ShadowMapping"] = _shadowMapping;
      _renderShader.attach();

      if (_analygraphMode)
	{
	  const double eyedist = 6.5; //
	  Vector eyeDisplacement(0.5 * eyedist, 0, 0);
	  
	  _viewPortInfo->buildMatrices(-eyeDisplacement);
	  _viewPortInfo->loadMatrices();
	  if (_shadowMapping)
	    _light0->loadShadowTextureMatrix(7);

	  glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
	  drawScene(*_renderTarget);
	  
	  _viewPortInfo->buildMatrices(eyeDisplacement);
	  _viewPortInfo->loadMatrices();
	  if (_shadowMapping)
	    _light0->loadShadowTextureMatrix(7);
	  
	  glClear(GL_DEPTH_BUFFER_BIT);
	  glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE);
	  drawScene(*_renderTarget);
	  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	}
      else
	{
	  _viewPortInfo->buildMatrices();
	  _viewPortInfo->loadMatrices();
	  if (_shadowMapping)
	    _light0->loadShadowTextureMatrix(7);
	  drawScene(*_renderTarget);
	}
      
      _renderTarget->detach();

      //////////////FILTERING////////////
      //Store what the last FBO was for later blitting to the screen
      magnet::GL::FBO* lastFBO = &(*_renderTarget);
      
      bool FBOalternate = false;

      if (_filterEnable)
	{
	  //Check if we need an extra pass where we calculate normals and depth values
	  bool renderNormsAndDepth = false;
	  
	  for (Gtk::TreeModel::iterator iPtr = _filterStore->children().begin(); 
	       iPtr != _filterStore->children().end(); ++iPtr)
	    {
	      void* filter_ptr = (*iPtr)[_filterModelColumns->m_filter_ptr];

	      if (static_cast<coil::filter*>(filter_ptr)->needsNormalDepth())
		{ renderNormsAndDepth = true; break; }
	    }

	  if (renderNormsAndDepth)
	    {
	      _normalsFBO.attach();
	      _nrmlShader.attach();
	      drawScene(_normalsFBO);
	      _normalsFBO.detach();
	    }

	  //Bind the original image to texture (unit 0)
	  _renderTarget->getColorTexture().bind(0);
	  
	  //Now bind the texture which has the normals and depths (unit 1)
	  _normalsFBO.getColorTexture().bind(1);

	  //High quality depth information is attached to (unit 2)
	  _renderTarget->getDepthTexture().bind(2);

	  for (Gtk::TreeModel::iterator iPtr = _filterStore->children().begin(); 
	       iPtr != _filterStore->children().end(); ++iPtr)
	    {
	      void* filter_ptr = (*iPtr)[_filterModelColumns->m_filter_ptr];
	      coil::filter& filter = *static_cast<coil::filter*>(filter_ptr);

	      if (!((*iPtr)[_filterModelColumns->m_active])) continue; //Only run active filters, skip to the next filter

	      if (filter.type_id() == coil::detail::filterEnum<coil::FlushToOriginal>::val)
		{//Check if we're trying to flush the drawing
		  lastFBO->attach();
		  glActiveTextureARB(GL_TEXTURE0);
		  //Now copy the texture 
		  glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, _viewPortInfo->getWidth(), _viewPortInfo->getHeight());
		  lastFBO->detach();
		}
	      else
		{
		  //The last output goes into texture 3
		  lastFBO->getColorTexture().bind(3);
		  
		  if (FBOalternate)
		    _filterTarget1.attach();
		  else
		    _filterTarget2.attach();
		  
		  filter.invoke(3, _viewPortInfo->getWidth(), _viewPortInfo->getHeight(), *_viewPortInfo);
		  
		  if (FBOalternate)
		    _filterTarget1.detach();
		  else
		    _filterTarget2.detach();
		  
		  lastFBO = FBOalternate ? &_filterTarget1 : &_filterTarget2;
		  
		  FBOalternate = !FBOalternate;
		}
	    }
	}

      //Restore the fixed pipeline
      //And turn off the shadow texture
      glUseProgramObjectARB(0);
      //Now blit the stored scene to the screen
      lastFBO->blitToScreen(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
    }
  else    
    {
      _renderTarget->attach();
      _viewPortInfo->loadMatrices();
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      drawScene(*_renderTarget);
      _renderTarget->detach();
      _renderTarget->blitToScreen(_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
    }
  
  //We clear the depth as merely disabling gives artifacts
  glClear(GL_DEPTH_BUFFER_BIT); 

  //Enter the interface draw for all objects
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->interfaceRender();

  glutSwapBuffers();

  //Check if we're recording and then check that if we're
  //framelocking, check that new data is available
  if ((_snapshot || _record) && (!_simframelock || _newData))
    {
      _newData = false;

      std::vector<magnet::image::Pixel<magnet::image::RGB> > pixels;
      pixels.resize(_viewPortInfo->getWidth() * _viewPortInfo->getHeight());
      //Read the pixels into our container
      glReadPixels(0,0, _viewPortInfo->getWidth(), _viewPortInfo->getHeight(), GL_RGB, 
		   GL_UNSIGNED_BYTE, &pixels[0]);
      
      std::string path;
      {
	Gtk::FileChooserButton* fileChooser;
	_refXml->get_widget("snapshotDirectory", fileChooser);
	path = fileChooser->get_filename();
      }

      if (_snapshot)
	{
	  _snapshot = false;

	  if (_PNGFileFormat)
	    magnet::image::writePNGFile(path + "/snapshot.png", pixels, _viewPortInfo->getWidth(), 
					_viewPortInfo->getHeight(), 9, false, true);
	  else
	    magnet::image::writeBMPFile(path + "/snapshot.bmp", pixels, _viewPortInfo->getWidth(), 
					_viewPortInfo->getHeight());
	}

      if (_record)
	{
	  std::ostringstream filename;
	  filename << std::setw(6) <<  std::setfill('0') << std::right << std::dec << _snapshot_counter++;
	  
	  if (_PNGFileFormat)
	    magnet::image::writePNGFile(path + "/" + filename.str() +".png", pixels, 
					_viewPortInfo->getWidth(), _viewPortInfo->getHeight(), 1, true, true);
	  else
	    magnet::image::writeBMPFile(path + "/" + filename.str() +".bmp", pixels, 
					_viewPortInfo->getWidth(), _viewPortInfo->getHeight());
	}                         
    }

  ++_frameCounter; 
  _lastFrameTime = _currFrameTime;
  _frameRenderTime = glutGet(GLUT_ELAPSED_TIME) - _currFrameTime;
}

void 
CLGLWindow::drawScene(magnet::GL::FBO& fbo)
{
  _light0->glUpdateLight();
  
  //Enter the render ticks for all objects
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->glRender(fbo);
  
  if (_showLight) _light0->drawLight();
}

void CLGLWindow::CallBackReshapeFunc(int w, int h)
{
  if (!CoilRegister::getCoilInstance().isRunning() || !_readyFlag) return;

  _viewPortInfo->setHeightWidth(h,w);
  //Setup the viewport
  glViewport(0, 0, w, h); 

  //Update the viewport
  _viewPortInfo->buildMatrices();
  _renderTarget->resize(w, h);
  
  if (_shaderPipeline)
    {
      _filterTarget1.resize(w, h);
      _filterTarget2.resize(w, h);
      _normalsFBO.resize(w, h);
    }

  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->resize(w, h);
  
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
	  if (_pickingEnabled) performPicking(x,y);
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
      _viewPortInfo->mouseMovement(diffX, diffY);
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

bool 
CLGLWindow::simupdateTick()
{
  ++_updateCounter;//For the updates per second

  for (;;)
    {
      _systemQueue->drainQueue();

      //Block the simulation if _simrun is false or if we're in frame lock
      //and a new frame has not been drawn.
      if (_simrun && (!_simframelock || (_lastUpdateTime != getLastFrameTime()))) break;
      
      //Jump out without an update if the window has been killed
      if (!isReady()) return false;

      //1ms delay to lower CPU usage while blocking, but not to affect framelocked render rates
      timespec sleeptime;
      sleeptime.tv_sec = 0;
      sleeptime.tv_nsec = 1000000;
      nanosleep(&sleeptime, NULL);
    }

  //Only redraw if the screen has actually refreshed
  if (_lastUpdateTime == getLastFrameTime()) return false;

  _lastUpdateTime = getLastFrameTime();

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

  _record = recordButton->get_active();  
}

void 
CLGLWindow::lightShowCallback()
{
  Gtk::CheckButton* lightShowButton;
  _refXml->get_widget("lightShow", lightShowButton);
  
  _showLight = lightShowButton->get_active();
}

void 
CLGLWindow::lightPlaceCallback()
{
  *_light0 = *_viewPortInfo;
}

void 
CLGLWindow::shadowIntensityCallback(double val)
{
  _shadowIntensity = val;
}

void 
CLGLWindow::snapshotFileFormatCallback()
{
  Gtk::RadioButton* radioButton;
  _refXml->get_widget("snapshotPNG", radioButton);
  _PNGFileFormat = radioButton->get_active();
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
  delete static_cast<coil::filter*>(tmp_ptr);

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
    [coil::filter::getSelectColumnsInstance().m_col_id];

  (*iter)[_filterModelColumns->m_filter_ptr] 
    = coil::filter::createFromID(type_id);

  (*iter)[_filterModelColumns->m_name]
    = coil::filter::getName(type_id);
  
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

      coil::filter* filter_ptr
	= (coil::filter*)((void*)((*iter)[_filterModelColumns->m_filter_ptr]));
      
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
      
      coil::filter* filter_ptr
	= (coil::filter*)((void*)((*iter)[_filterModelColumns->m_filter_ptr]));
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
	    delete static_cast<coil::filter*>(tmp_ptr);
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

extern const guint8 coilsplash[];
extern const size_t coilsplash_size;

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
CLGLWindow::renderNormalsCallback()
{
  Gtk::CheckButton* Button;
  _refXml->get_widget("showNormals", Button);
  
  bool showNormals = Button->get_active();

  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->setDisplayNormals(showNormals);
}

void
CLGLWindow::performPicking(int x, int y)
{
  glUseProgramObjectARB(0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);       
  _viewPortInfo->loadMatrices();
  //Perform unique coloring of screen objects

  cl_uint startVal = 0;
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->initPicking(startVal);

  //Now render the scene
  glDisable(GL_BLEND); 
  glDisable(GL_DITHER); 
  glDisable(GL_FOG); 
  glDisable(GL_LIGHTING); 
  glDisable(GL_TEXTURE_1D); 
  glDisable(GL_TEXTURE_2D); 
  glDisable(GL_TEXTURE_3D); 
  glShadeModel (GL_FLAT);

  //Enter the render ticks for all objects
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->pickingRender();

  unsigned char pixel[4];
  
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel);

  glEnable(GL_BLEND); 
  glEnable(GL_DITHER); 
  glEnable(GL_FOG); 
  glEnable(GL_LIGHTING); 
  glEnable(GL_TEXTURE_1D); 
  glEnable(GL_TEXTURE_2D); 
  glEnable(GL_TEXTURE_3D); 
  glShadeModel(GL_SMOOTH);

  //Now let the objects know what was picked
  const cl_uint objID = pixel[0] + 256 * (pixel[1] + 256 * (pixel[2] + 256 * pixel[3]));
  startVal = 0;
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    (*iPtr)->finishPicking(startVal, objID);
}

void 
CLGLWindow::rebuildRenderView()
{
  _renderObjStore->clear();
  
  for (std::vector<magnet::thread::RefPtr<RenderObj> >::iterator iPtr = RenderObjects.begin();
       iPtr != RenderObjects.end(); ++iPtr)
    {
      Gtk::TreeModel::iterator iter = _renderObjStore->append();
      
      (*iter)[_renderObjModelColumns->m_name]
	= (*iPtr)->getName();
      
      (*iter)[_renderObjModelColumns->m_visible]
	= (*iPtr)->isVisible();

      (*iter)[_renderObjModelColumns->m_id]
	= iPtr - RenderObjects.begin();
    }
}

void CLGLWindow::visibleRObjCallback() 
{
  Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
    _renderObjView->get_selection();
  Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

  Gtk::ToggleButton *visibleBtn;
  _refXml->get_widget("robjVisible", visibleBtn);

  if(iter)
    {
      size_t objID = (*iter)[_renderObjModelColumns->m_id]; 
      bool newState = visibleBtn->get_active();

      RenderObjects[objID]->setVisible(newState);
      (*iter)[_renderObjModelColumns->m_visible] = newState;
    }

  selectRObjCallback();
}
void CLGLWindow::editRObjCallback() {}
void CLGLWindow::deleteRObjCallback() {}
void CLGLWindow::addRObjCallback() {}
void CLGLWindow::selectRObjCallback() 
{
  Glib::RefPtr<Gtk::TreeSelection> refTreeSelection =
    _renderObjView->get_selection();

  Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();

  Gtk::Button *deleteBtn, *editBtn, *addBtn;
  Gtk::ToggleButton *visibleBtn;
  Gtk::Image *visibleImg;
  Gtk::ScrolledWindow* win;
  _refXml->get_widget("robjDelete", deleteBtn);
  _refXml->get_widget("robjEdit", editBtn);
  _refXml->get_widget("robjAdd", addBtn);
  _refXml->get_widget("robjVisible", visibleBtn);
  _refXml->get_widget("robjVisibleImg", visibleImg);
  _refXml->get_widget("ObjectOptions", win);

  win->remove(); //Clear the current object controls
  if(iter)
    {
      //Enable the filter buttons
      deleteBtn->set_sensitive(false);
      editBtn->set_sensitive(false); 
      visibleBtn->set_sensitive(true);

      if (RenderObjects[(*iter)[_renderObjModelColumns->m_id]]->isVisible())
	{//Object is visible
	  visibleBtn->set_active(true);
	  visibleImg->set(Gtk::Stock::YES, Gtk::ICON_SIZE_BUTTON);
	}
      else
	{//Object is not visible
	  visibleBtn->set_active(false);
	  visibleImg->set(Gtk::Stock::NO, Gtk::ICON_SIZE_BUTTON);
	}

      //Load the controls for the window
      RenderObjects[(*iter)[_renderObjModelColumns->m_id]]->showControls(win);
    }
  else
    {
      //Disable all of the filter buttons
      deleteBtn->set_sensitive(false);
      editBtn->set_sensitive(false); 
      visibleBtn->set_sensitive(false);
    }

  addBtn->set_sensitive(false); 
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
    _light0->setFOVY(FOVscale->get_value());
    _light0->buildMatrices();
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
    _viewPortInfo->setSimUnitLength(boost::lexical_cast<double>(val));
  }

  {
    Gtk::Entry* pixelPitch;
    _refXml->get_widget("pixelPitch", pixelPitch);
    std::string val = pixelPitch->get_text();
    if (val.empty()) {val = "0.25"; pixelPitch->set_text("0.25"); }
    _viewPortInfo->setPixelPitch(boost::lexical_cast<double>(val) / 10);
  }

  {
    Gtk::Label* XHead;
    _refXml->get_widget("XHead", XHead);
    Gtk::Label* YHead;
    _refXml->get_widget("YHead", YHead);
    Gtk::Label* ZHead;
    _refXml->get_widget("ZHead", ZHead);
    std::ostringstream os;
    os << _viewPortInfo->getHeadLocation()[0] << "cm";
    XHead->set_text(os.str());
    os.str("");
    os << _viewPortInfo->getHeadLocation()[1] << "cm";
    YHead->set_text(os.str());
    os.str("");
    os << _viewPortInfo->getHeadLocation()[2] << "cm";
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
  _viewPortInfo->setHeadLocation(Vector(0,0,_viewPortInfo->getHeadLocation()[2]));
  _viewPortInfo->setFOVY(60.f, false);
}

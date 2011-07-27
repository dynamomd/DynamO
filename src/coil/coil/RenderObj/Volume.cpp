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

#include "Volume.hpp"
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <magnet/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

namespace coil {
  void 
  RVolume::deinit()
  {
    _currentDepthFBO.deinit();
    _data.deinit();
    _transferFuncTexture.deinit();
    _shader.deinit();
    _cube.deinit();
  }

  void 
  RVolume::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue) 
  {
    RenderObj::init(systemQueue);
    _shader.build();
    _cube.init();
    _transferFuncTexture.init(256);
    initGTK();
  }

  void 
  RVolume::loadRawFile(std::string filename, size_t width, size_t height, size_t depth,
		       size_t bytes)
  {
    std::ifstream file(filename.c_str(), std::ifstream::binary);
    std::vector<GLubyte> inbuffer(width * height * depth);
    
    switch (bytes)
      {
      case 1:
	{
	  file.read(reinterpret_cast<char*>(&inbuffer[0]), inbuffer.size());
	  if (file.fail()) M_throw() << "Failed to load the texture from the file";
	}
	break;
      case 2:
	{
	  std::vector<uint16_t> tempBuffer(width * height * depth);
	  file.read(reinterpret_cast<char*>(&tempBuffer[0]), 2 * tempBuffer.size());
	  if (file.fail()) M_throw() << "Failed to load the texture from the file";
	  for (size_t i(0); i < tempBuffer.size(); ++i)
	    inbuffer[i] = uint8_t(tempBuffer[i] >> 8);
	}
	break;
      default:
	M_throw() << "Cannot load at that bit depth yet";
      }

    loadData(inbuffer, width, height, depth);
  }

  void
  RVolume::loadSphereTestPattern()
  {
    const size_t size(256);

    std::vector<GLubyte> inbuffer(size * size * size);

    //Sphere test pattern
    for (size_t z(0); z < size; ++z)
      for (size_t y(0); y < size; ++y)
        for (size_t x(0); x < size; ++x)
          inbuffer[x + size * (y + size * z)] 
	    = (std::sqrt(std::pow(x - size / 2.0, 2) 
			 + std::pow(y - size / 2.0, 2) 
			 + std::pow(z - size / 2.0, 2))
	       < 122.0) ? 255.0 : 0;
    
    loadData(inbuffer, size, size, size);
  }

  namespace {
    inline size_t coordCalc(GLint x, GLint y, GLint z, 
			    GLint width, GLint height, GLint depth)
    {
      x = magnet::clamp(x, 0, width  - 1);
      y = magnet::clamp(y, 0, height - 1);
      z = magnet::clamp(z, 0, depth  - 1);
      return x + width * (y + height * z);
    }
  }

  void 
  RVolume::loadData(const std::vector<GLubyte>& inbuffer, size_t width, size_t height, size_t depth)
  {
    std::vector<GLubyte> voldata(4 * width * height * depth);
    
    std::vector<float>& histogram = _transferFunction->getHistogram();
    histogram = std::vector<float>(256, 0);
    
    for (int z(0); z < int(depth); ++z)
      for (int y(0); y < int(height); ++y)
	for (int x(0); x < int(width); ++x)
	  {
	    Vector sample1(inbuffer[coordCalc(x - 1, y, z, width, height, depth)],
			   inbuffer[coordCalc(x, y - 1, z, width, height, depth)],
			   inbuffer[coordCalc(x, y, z - 1, width, height, depth)]);
	    
	    Vector sample2(inbuffer[coordCalc(x + 1, y, z, width, height, depth)],
			   inbuffer[coordCalc(x, y + 1, z, width, height, depth)],
			   inbuffer[coordCalc(x, y, z + 1, width, height, depth)]);
	    
	    //Note, we store the negative gradient (we point down
	    //the slope)
	    Vector grad = sample1 - sample2;
	    
	    float nrm = grad.nrm();
	    if (nrm > 0) grad /= nrm;
	    
	    size_t coord = x + width * (y + height * z);
	    voldata[4 * coord + 0] = uint8_t((grad[0] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 1] = uint8_t((grad[1] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 2] = uint8_t((grad[2] * 0.5 + 0.5) * 255);
	    
	    GLubyte val = inbuffer[coordCalc(x, y, z, width, height, depth)];
	    voldata[4 * coord + 3] 
	      = val;

	    histogram[val] += 1;
	  }
    
    {
      float logMaxVal = std::log(*std::max_element(histogram.begin(), histogram.end()));
      float logMinVal = std::log(std::max(*std::min_element(histogram.begin(), histogram.end()), 1.0f));
      float normalization = 1.0 / (logMaxVal - logMinVal);

      for (std::vector<float>::iterator iPtr = histogram.begin();
	   iPtr != histogram.end(); ++iPtr)
	{
	  if (*iPtr == 0) *iPtr = 1.0;
	  *iPtr = (std::log(*iPtr) - logMinVal) * normalization;
	}
    }

    _data.init(width, height, depth);
    _data.subImage(voldata, GL_RGBA);
  }
  
  void 
  RVolume::glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& camera)
  {
    if (!_visible || !_data.isValid()) return;

    //Resize the copy FBO
    _currentDepthFBO.init(fbo.getWidth(), fbo.getHeight());
    //Before we render, we need the current depth buffer for depth testing.
    fbo.detach();
    fbo.copyto(_currentDepthFBO, GL_DEPTH_BUFFER_BIT);
    fbo.attach();

    //Now bind this copied depth texture to texture unit 0
    _currentDepthFBO.getDepthTexture().bind(0);
    _data.bind(1);
    _transferFuncTexture.bind(2);

    //You must save the current shader before binding the volume
    //shader
    GLhandleARB oldshader = glGetHandleARB(GL_PROGRAM_OBJECT_ARB);


    _shader["ProjectionMatrix"] = camera.getProjectionMatrix();
    _shader["ViewMatrix"] = camera.getViewMatrix();
    _shader["FocalLength"] = GLfloat(1.0f / std::tan(camera.getFOVY() * (M_PI / 360.0f)));
    { 
      std::tr1::array<GLfloat,2> winsize = {{camera.getWidth(), camera.getHeight()}};
      _shader["WindowSize"] = winsize;
    }
    _shader["RayOrigin"] = camera.getEyeLocation();
    _shader["DepthTexture"] = 0;
    _shader["NearDist"] = camera.getZNear();
    _shader["FarDist"] = camera.getZFar();
    _shader["DataTexture"] = 1;
    _shader["StepSize"] = _stepSizeVal;
    _shader["DiffusiveLighting"] = GLfloat(_diffusiveLighting->get_value());
    _shader["SpecularLighting"] = GLfloat(_specularLighting->get_value());
    _shader["DitherRay"] = GLfloat(_ditherRay->get_value());
    _shader["TransferTexture"] = 2;
    _shader["LightPosition"] = Vector(-2,-2,-2);
    _shader.attach();
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);

    _currentDepthFBO.getContext().color(1,1,0,1);
    _cube.glRender();

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);

    glUseProgramObjectARB(oldshader);
  }

  void 
  RVolume::transferFunctionUpdated()
  {
    if (_transferFunction.get() != NULL)
      _transferFuncTexture.subImage(_transferFunction->getColorMap(), GL_RGBA);
  }
  
  void
  RVolume::initGTK()
  {
    _optList.reset(new Gtk::VBox);//The Vbox of options   

    {//Transfer function widget
      _transferFunction.reset(new magnet::gtk::TransferFunction
			      (magnet::function::MakeDelegate
			       (this, &RVolume::transferFunctionUpdated)));
      _transferFunction->set_size_request(-1, 100);
      
      _optList->add(*_transferFunction); _transferFunction->show();
      transferFunctionUpdated(); //Force an update of the transfer function now we have the widget
    }

    {//Volume renderer step size
      _stepSize.reset(new Gtk::Entry);
      Gtk::HBox* box = manage(new Gtk::HBox);	
      Gtk::Label* label = manage(new Gtk::Label("Raytrace Step Size"));
      box->pack_start(*label, false, false); label->show();
      box->pack_end(*_stepSize, false, false);
      _stepSize->show(); _stepSize->set_text("0.01");      
      _optList->add(*box); box->show();
    }

    {//Diffusive lighting
      Gtk::HBox* box = manage(new Gtk::HBox);	
      Gtk::Label* label = manage(new Gtk::Label("Diffuse Lighting"));
      box->pack_start(*label, false, false); label->show();
      _diffusiveLighting.reset(new Gtk::HScale);
      box->pack_end(*_diffusiveLighting, true, true);
      _diffusiveLighting->set_range(0,2);
      _diffusiveLighting->set_digits(3);
      _diffusiveLighting->show();      
      _diffusiveLighting->set_value(1.0);
      _optList->add(*box); box->show();
    }

    {//Specular lighting
      Gtk::HBox* box = manage(new Gtk::HBox);	
      Gtk::Label* label = manage(new Gtk::Label("Specular Lighting"));
      box->pack_start(*label, false, false); label->show();
      _specularLighting.reset(new Gtk::HScale);
      box->pack_end(*_specularLighting, true, true);
      _specularLighting->set_range(0,2);
      _specularLighting->set_digits(3);
      _specularLighting->show();      
      _specularLighting->set_value(1.0);
      _optList->add(*box); box->show();
    }

    {//Ray Dithering
      Gtk::HBox* box = manage(new Gtk::HBox);	
      Gtk::Label* label = manage(new Gtk::Label("Ray Dithering"));
      box->pack_start(*label, false, false); label->show();
      _ditherRay.reset(new Gtk::HScale);
      box->pack_end(*_ditherRay, true, true);
      _ditherRay->set_range(0, 1);
      _ditherRay->set_digits(3);
      _ditherRay->show();
      _ditherRay->set_value(1.0);
      _optList->add(*box); box->show();
    }
    
    _optList->show();
    //Callbacks
    _stepSize->signal_changed()
      .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_stepSize));
    _stepSize->signal_activate().connect(sigc::mem_fun(*this, &RVolume::guiUpdate));

    guiUpdate();
  }

  void
  RVolume::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _optList->unparent();
    win->add(*_optList);
    win->show();
  }

  void 
  RVolume::guiUpdate()
  {
    std::string val = _stepSize->get_text();
    if (val.empty()) {val = "0.01"; _stepSize->set_text("0.01"); }
    
    _stepSizeVal = boost::lexical_cast<double>(val);
  }
}

//////////////////SOME OLD CODE THAT MIGHT BE USEFUL
//      //Gaussian blur the system
//      const double weights[5] = {0.0544886845,0.244201342,0.4026199469,0.244201342,0.0544886845};
//
//      for (size_t component(0); component < 3; ++component)
//	for (GLint z(0); z < _depth; ++z)
//	  for (GLint y(0); y < _height; ++y)
//	    for (GLint x(0); x < _width; ++x)
//	      {
//		double sum(0);
//		for (int i(0); i < 5; ++i)
//		  {
//		    GLint pos[3] = {x,y,z};
//		    pos[component] += i - 2;
//		    sum += weights[i] 
//		      * voldata[4 * detail::coordCalc(pos[0], pos[1], pos[2], 
//						      _width, _height, _depth) + component];
//		  }
//		voldata[4 * detail::coordCalc(x, y, z, _width, _height, _depth) 
//			+ component] = sum;
//	      }
//      
//      //Renormalize the gradients
//      for (GLint z(0); z < _depth; ++z)
//	for (GLint y(0); y < _height; ++y)
//	  for (GLint x(0); x < _width; ++x)
//	    {
//	      std::vector<GLubyte>::iterator iPtr = voldata.begin()
//		+ 4 * detail::coordCalc(x, y, z, _width, _height, _depth);
//	      
//	      Vector grad(*(iPtr + 0) - 128.0, *(iPtr + 1) - 128.0, *(iPtr + 2) - 128.0);
//	      float nrm = grad.nrm();
//	      if (nrm > 0) grad /= nrm;
//
//	      *(iPtr + 0) = uint8_t((grad[0] * 0.5 + 0.5) * 255);
//	      *(iPtr + 1) = uint8_t((grad[1] * 0.5 + 0.5) * 255);
//	      *(iPtr + 2) = uint8_t((grad[2] * 0.5 + 0.5) * 255);
//	    }

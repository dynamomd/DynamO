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

#include <coil/RenderObj/Volume.hpp>
#include <coil/RenderObj/Light.hpp>
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <magnet/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

extern const guint8 Volume_Icon[];
extern const size_t Volume_Icon_size;

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RVolume::getIcon()
  {
    return Gdk::Pixbuf::create_from_inline(Volume_Icon_size, Volume_Icon);
  }

  void 
  RVolume::deinit()
  {
    _currentDepthFBO.deinit();
    _data.deinit();
    _transferFuncTexture.deinit();
    _shader.deinit();
    _depthCopyShader.deinit();
    _cube.deinit();
  }

  void 
  RVolume::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue) 
  {
    RenderObj::init(systemQueue);
    _shader.build();
    _depthCopyShader.build();
    _cube.init();
    _transferFuncTexture.init(256, GL_RGBA16F);
    _transferFuncTexture.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    _transferFuncTexture.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    _transferFuncTexture.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

    _preintTransferFuncTexture.init(256, GL_RGBA16F);
    _preintTransferFuncTexture.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    _preintTransferFuncTexture.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    _preintTransferFuncTexture.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

    //Resize the copy FBO
    //Build depth buffer
    std::tr1::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
    depthTexture->init(800, 600, GL_DEPTH_COMPONENT);
    depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
    _currentDepthFBO.init();
    _currentDepthFBO.attachTexture(depthTexture);

    initGTK();
  }

  void 
  RVolume::loadRawFile(std::string filename, size_t width, size_t height, 
		       size_t depth, size_t bytes)
  {
    std::tr1::array<size_t, 3> dim = {{width, height, depth}};
    _currentDepthFBO.getContext().queueTask(magnet::function::Task::makeTask
					    (&RVolume::loadRawFileWorker, this, 
					     filename, dim, bytes));
  }

  void 
  RVolume::loadRawFileWorker(std::string filename, std::tr1::array<size_t,3> dim, 
			     size_t bytes)
  {
    std::ifstream file(filename.c_str(), std::ifstream::binary);
    std::vector<GLubyte> inbuffer(dim[0] * dim[1] * dim[2]);
    
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
	  std::vector<uint16_t> tempBuffer(dim[0] * dim[1] * dim[2]);
	  file.read(reinterpret_cast<char*>(&tempBuffer[0]), 2 * tempBuffer.size());
	  if (file.fail()) M_throw() << "Failed to load the texture from the file";
	  for (size_t i(0); i < tempBuffer.size(); ++i)
	    inbuffer[i] = uint8_t(tempBuffer[i] >> 8);
	}
	break;
      default:
	M_throw() << "Cannot load at that bit depth yet";
      }

    //Debug loading of data
    //loadSphereTestPattern();

    loadData(inbuffer, dim[0], dim[1], dim[2]);
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
	    = std::sqrt(std::pow(x - size / 2.0, 2)
			 + std::pow(y - size / 2.0, 2) 
			 + std::pow(z - size / 2.0, 2))
	    ;
    
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
	    
	    //Do a central difference scheme
	    Vector grad = sample1 - sample2;
	    
	    float nrm = grad.nrm();
	    if (nrm > 0) grad /= nrm;

	    size_t coord = x + width * (y + height * z);
	    voldata[4 * coord + 0] = uint8_t((grad[0] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 1] = uint8_t((grad[1] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 2] = uint8_t((grad[2] * 0.5 + 0.5) * 255);
	    
	    GLubyte val = inbuffer[coordCalc(x, y, z, width, height, depth)];
	    voldata[4 * coord + 3] = val;
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
    _data.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    _data.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    _data.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    _data.parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    _data.parameter(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    _data.subImage(voldata, GL_RGBA);
  }
  
  void 
  RVolume::forwardRender(magnet::GL::FBO& fbo,
			 const magnet::GL::Camera& camera,
			 const RLight& light,
			 RenderMode mode)
  {
    if (!_visible || !_data.isValid()) return;

    //Before we render, we need the current depth buffer so we can test against it
//    fbo.detach();
//    _currentDepthFBO.resize(fbo.getWidth(), fbo.getHeight());
//    _currentDepthFBO.attach();
//    glClear(GL_DEPTH_BUFFER_BIT);
//    _depthCopyShader.attach();
//    fbo.getDepthTexture()->bind(0);
//    _depthCopyShader["depthTex"] = 0;
//    _depthCopyShader.invoke();
//    _depthCopyShader.detach();
//    _currentDepthFBO.detach();
//    fbo.attach();

    //Now bind this copied depth texture to texture unit 0
    _currentDepthFBO.getDepthTexture()->bind(0);
    _data.bind(1);
    _transferFuncTexture.bind(2);
    _preintTransferFuncTexture.bind(3);

    _shader.attach();
    _shader["ProjectionMatrix"] = camera.getProjectionMatrix();
    _shader["ViewMatrix"] = camera.getViewMatrix();
    _shader["FocalLength"] = GLfloat(1.0 / std::tan(camera.getFOVY() * (M_PI / 360.0)));
    { 
      std::tr1::array<GLfloat,2> winsize = {{camera.getWidth(), camera.getHeight()}};
      _shader["WindowSize"] = winsize;
    }
    _shader["RayOrigin"] = camera.getEyeLocationObjSpace();
    _shader["DepthTexture"] = 0;
    _shader["DataTexture"] = 1;
    _shader["StepSize"] = _stepSizeVal;

    _shader["lightColor"] = light.getColor();
    _shader["lightAttenuation"] = light.getAttenuation();
    _shader["lightSpecularExponent"] = light.getSpecularExponent();
    _shader["lightSpecularFactor"] = light.getSpecularFactor();
    _shader["lightIntensity"] = light.getIntensity();
    _shader["lightPosition"] = light.getEyespacePosition(camera);

    _shader["DitherRay"] = GLint(_ditherRay->get_active());
    _shader["TransferTexture"] = 2;
    _shader["IntTransferTexture"] = 3;
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);

    _currentDepthFBO.getContext().cleanupAttributeArrays();
    _cube.glRender();
    _shader.detach();

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
  }

  void 
  RVolume::transferFunctionUpdated()
  {
    if (_transferFunction.get() != NULL)
      {
	size_t samples = 256;
	float transmittanceFactor = 1000;

	std::vector<float> data = _transferFunction->getMap(samples, transmittanceFactor);
	std::vector<GLfloat> GLdata = data;
	_transferFuncTexture.subImage(GLdata, GL_RGBA);
	
	data = _transferFunction->getPreIntegratedMap(samples, transmittanceFactor);
	GLdata = data;
	_preintTransferFuncTexture.subImage(GLdata, GL_RGBA);
      }
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

    {//Ray Dithering and filtering
      Gtk::HBox* box = manage(new Gtk::HBox);
      _ditherRay.reset(new Gtk::CheckButton("Dither"));
      _filterData.reset(new Gtk::CheckButton("Filter Data"));
      
      _ditherRay->set_active(true);
      _ditherRay->show();
      _filterData->set_active(true);
      _filterData->show();

      box->pack_end(*_ditherRay, true, true);
      box->pack_end(*_filterData, true, true);
      _optList->add(*box); box->show();
    }
    
    _optList->show();
    //Callbacks
    _stepSize->signal_changed()
      .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_stepSize));
    _stepSize->signal_activate().connect(sigc::mem_fun(*this, &RVolume::guiUpdate));

    _filterData->signal_toggled()
      .connect(sigc::mem_fun(*this, &RVolume::guiUpdate));

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

    if (_filterData->get_active())
      {
	_data.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	_data.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }
    else
      {
	_data.parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	_data.parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }    
  }
}

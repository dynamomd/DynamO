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

#include <coil/RenderObj/Volume.hpp>
#include <coil/RenderObj/Light.hpp>
#include <coil/RenderObj/console.hpp>
#include <magnet/GL/objects/primitives/cube.hpp>

#ifdef COIL_TIFFSUPPORT
# include <magnet/image/TIFF.hpp>
#endif

#include <magnet/gtk/numericEntry.hpp>
#include <magnet/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

#include <coil/images/images.hpp>

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RVolume::getIcon()
  { return coil::images::Volume_Icon(); }

  void 
  RVolume::deinit()
  {
    _currentDepthFBO.deinit();
    _data.deinit();
    _transferFuncTexture.deinit();
    _shader.deinit();
    _depthCopyShader.deinit();
    _cubeVertices.deinit();
  }

  void 
  RVolume::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue);
    _shader.defines("LIGHT_COUNT") = 1;
    _shader.build();
    _depthCopyShader.build();
    _cubeVertices.init(magnet::GL::objects::primitives::Cube::getVertices(), 3);
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
    std::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D);
    depthTexture->init(800, 600, GL_DEPTH_COMPONENT);
    depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
    _currentDepthFBO.init();
    _currentDepthFBO.attachTexture(depthTexture);

    initGTK();
    _initialised = true;
  }

  void 
  RVolume::loadRawFile(std::string filename, std::array<size_t, 3> dim, size_t bytes)
  {
    std::ifstream file(filename.c_str(), std::ifstream::binary);
    std::vector<uint8_t> filebuffer(dim[0] * dim[1] * dim[2] * bytes);
    file.read(reinterpret_cast<char*>(&filebuffer[0]), filebuffer.size());
    if (file.fail()) M_throw() << "Failed to load the texture from the file, possible incorrect dimensions ";
	  
    //Debug loading of data
    //loadSphereTestPattern();

    std::vector<GLubyte> outbuffer(dim[0] * dim[1] * dim[2]);
    for (size_t x(0); x < dim[0]; ++x)
      for (size_t y(0); y < dim[1]; ++y)
	for (size_t z(0); z < dim[2]; ++z)
	  outbuffer[x + (y + z* dim[1]) * dim[0]] = filebuffer[(x + (y + z * dim[1]) * dim[0]) * bytes];

    //size_t maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    //_dimensions = Vector{double(dim[0]) / maxdim, double(dim[1]) / maxdim, double(dim[2]) / maxdim};

    loadData(outbuffer, dim, Vector{1,1,1});
  }

#ifdef COIL_TIFFSUPPORT
  void 
  RVolume::loadTIFFFiles(std::vector<std::string> files)
  {
    const auto data = magnet::image::loadTIFFStack(files);

    std::array<size_t, 3> dim;
    dim[0] = data.width;
    dim[1] = data.height;
    dim[2] = data.depth;
    std::vector<GLubyte> outbuffer(data.width * data.height * data.depth);
    for (size_t x(0); x < data.width; ++x)
      for (size_t y(0); y < data.height; ++y)
	for (size_t z(0); z < data.depth; ++z)
	  {
	    const size_t idx = x + data.width * (y + data.height * z);
	    outbuffer[idx] = data.pixels[idx].r;
	  }

    loadData(outbuffer, dim, Vector{1,1,1});
  }
#endif

  void
  RVolume::loadSphereTestPattern()
  {
    const size_t size(256);

    std::vector<GLubyte> inbuffer(size * size * size);

    //Sphere test pattern
    for (size_t z(0); z < size; ++z)
      for (size_t y(0); y < size; ++y)
        for (size_t x(0); x < size; ++x)
          inbuffer[x + size * (y + size * z)] = std::sqrt(std::pow(x - size / 2.0, 2) + std::pow(y - size / 2.0, 2) + std::pow(z - size / 2.0, 2));
    
    loadData(inbuffer, std::array<size_t, 3>{{size, size, size}}, Vector{1,1,1});
  }

  namespace {
    inline GLubyte coordCalc(GLint x, GLint y, GLint z, 
			     GLint width, GLint height, GLint depth, const std::vector<GLubyte>& buffer)
    {
      if (((x < 0) || (x >= width)) 
	  || ((y < 0) || (y >= height))
	  || ((z < 0) || (z >= depth)))
	return 0;
      else
	return buffer[x + width * (y + height * z)];
    }
  }

  
  void
  RVolume::loadData(const std::vector<GLubyte>& inbuffer, std::array<size_t, 3> dim, Vector dimensions)
  {
    _dimensions = dimensions;
    
    //Figure out what the minimum step size is to capture all detail
    //of the model (nyquist sampling rate)
    _stepSizeVal = HUGE_VAL;
    for (size_t i(0); i < 3; ++i)
      _stepSizeVal = std::min(_stepSizeVal, GLfloat(0.5 * dimensions[i] / dim[i]));
    if (_stepSize)
      _stepSize->set_text(boost::lexical_cast<std::string>(_stepSizeVal));

    std::vector<GLubyte> voldata(4 * dim[0] * dim[1] * dim[2]);
    
    std::vector<float>& histogram = _transferFunction->getHistogram();
    histogram = std::vector<float>(256, 0);

    size_t width = dim[0];
    size_t height = dim[1];
    size_t depth = dim[2];

    for (int z(0); z < int(depth); ++z)
      for (int y(0); y < int(height); ++y)
	for (int x(0); x < int(width); ++x)
	  {
	    Vector sample1{double(coordCalc(x - 1, y, z, width, height, depth, inbuffer)),
		double(coordCalc(x, y - 1, z, width, height, depth, inbuffer)),
		double(coordCalc(x, y, z - 1, width, height, depth, inbuffer))};
	    
	    Vector sample2{double(coordCalc(x + 1, y, z, width, height, depth, inbuffer)),
		double(coordCalc(x, y + 1, z, width, height, depth, inbuffer)),
		double(coordCalc(x, y, z + 1, width, height, depth, inbuffer))};
	    
	    //Do a central difference scheme
	    Vector grad = sample1 - sample2;
	    
	    float nrm = grad.nrm();
	    if (nrm > 0) grad /= nrm;

	    size_t coord = x + width * (y + height * z);
	    voldata[4 * coord + 0] = uint8_t((grad[0] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 1] = uint8_t((grad[1] * 0.5 + 0.5) * 255);
	    voldata[4 * coord + 2] = uint8_t((grad[2] * 0.5 + 0.5) * 255);
	    
	    GLubyte val = coordCalc(x, y, z, width, height, depth, inbuffer);
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
			 std::vector<std::shared_ptr<RLight> >& lights,
			 GLfloat ambient,
			 RenderMode mode)
  {
    if (!_visible || !_data.isValid()) return;
    if (lights.empty()) return;

    //Before we render, we need the current depth buffer so we can test against it
    fbo.detach();

    if ((fbo.getWidth() != _currentDepthFBO.getWidth())
	|| (fbo.getHeight() != _currentDepthFBO.getHeight()))
      {
	_currentDepthFBO.deinit();
	std::shared_ptr<magnet::GL::Texture2D> 
	  depthTexture(new magnet::GL::Texture2D);
	depthTexture->init(fbo.getWidth(), fbo.getHeight(), GL_DEPTH_COMPONENT);
	depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	//depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
	_currentDepthFBO.init();
	_currentDepthFBO.attachTexture(depthTexture);
      }

    _currentDepthFBO.attach();
    glClear(GL_DEPTH_BUFFER_BIT);
    _depthCopyShader.attach();
    fbo.getDepthTexture()->bind(0);
    _depthCopyShader["depthTex"] = 0;
    _depthCopyShader.invoke();
    _depthCopyShader.detach();
    _currentDepthFBO.detach();
    fbo.attach();

    //Now bind this copied depth texture to texture unit 0
    _currentDepthFBO.getDepthTexture()->bind(0);
    _data.bind(1);
    _transferFuncTexture.bind(2);
    _preintTransferFuncTexture.bind(3);

    _shader.defines("LIGHT_COUNT") = lights.size();

    _shader.attach();

    std::vector<Vector> light_positions;
    std::vector<Vector> light_color;
    std::vector<Vector> light_factors;
    for (std::vector<std::shared_ptr<RLight> >::const_iterator 
	   iPtr = lights.begin(); iPtr != lights.end(); ++iPtr)
      {
	light_positions.push_back((*iPtr)->getEyespacePosition(camera));
	auto color = (*iPtr)->getLightColor();
	Vector vcol{color[0], color[1], color[2]};
	light_color.push_back(vcol);
	light_factors.push_back(Vector{0.0, (*iPtr)->getSpecularExponent(),(*iPtr)->getSpecularFactor()});
      }

    _shader["lightPosition"] = light_positions;
    _shader["lightColor"] = light_color;
    _shader["lightFactors"] = light_factors;
    _shader["RayOrigin"] = camera.getPosition();
    _shader["TransferTexture"] = 2;
    _shader["IntTransferTexture"] = 3;
    _shader["DepthTexture"] = 0;
    _shader["DataTexture"] = 1;
    _shader["StepSize"] = _stepSizeVal;
    _shader["DitherRay"] = GLint(_ditherRay->get_active());
    _shader["ProjectionMatrix"] = camera.getProjectionMatrix();
    _shader["ViewMatrix"] = camera.getViewMatrix();

    Vector volumeMin = double(-0.5) * _dimensions;
    Vector volumeMax = double(+0.5) * _dimensions;

    Vector invVolumeDimensions = Vector{1 / (volumeMax[0] - volumeMin[0]),
					1 / (volumeMax[1] - volumeMin[1]),
					1 / (volumeMax[2] - volumeMin[2])};

    _shader["volumeMin"] = volumeMin;
    _shader["volumeMax"] = volumeMax;
    _shader["invVolumeDimensions"] = invVolumeDimensions;
    _shader["ambientLight"] = ambient;
    
    _currentDepthFBO.getContext().setCullFace(true);
    _currentDepthFBO.getContext().setDepthTest(false);
    glCullFace(GL_FRONT);
    glDepthMask(GL_FALSE);


    _currentDepthFBO.getContext().cleanupAttributeArrays();
    _currentDepthFBO.getContext().setAttribute(magnet::GL::Context::instanceScaleAttrIndex, 
					       _dimensions[0],
					       _dimensions[1],
					       _dimensions[2], 1);

    _cubeVertices.drawArray(magnet::GL::element_type::TRIANGLES);
    _shader.detach();

    _currentDepthFBO.getContext().setDepthTest(true);
    _currentDepthFBO.getContext().setCullFace(false);
    glDepthMask(GL_TRUE);
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
      _transferFunction.reset(new magnet::gtk::TransferFunction(magnet::Delegate<void()>::create<RVolume, &RVolume::transferFunctionUpdated>(this)));
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
      _stepSize->show(); 
      _stepSize->set_text(boost::lexical_cast<std::string>(_stepSizeVal));
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
      .connect(sigc::bind(&magnet::gtk::forceNumericEntry, _stepSize.get()));
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

    if (_data.isValid())
      {
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
}

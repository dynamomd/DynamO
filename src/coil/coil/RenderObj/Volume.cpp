/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include <magnet/color/transferFunction.hpp>
#include <magnet/PNG.hpp>
#include <magnet/clamp.hpp>
#include <boost/lexical_cast.hpp>

namespace coil {
  RVolume::RVolume(std::string name):
    RQuads(name),
    _stepSizeVal(0.01)
  {}

  RVolume::~RVolume()
  {}

  void 
  RVolume::releaseCLGLResources()
  {
    if (_fbo.get() != NULL) _fbo->deinit();
    _fbo.release();
    _data.deinit();
    _transferFuncTexture.deinit();
    _shader.deinit();
    RQuads::releaseCLGLResources();
  }

  void 
  RVolume::initOpenGL() 
  {
    _shader.build();
    _fbo.reset(new magnet::GL::FBO);
    _fbo->init(_viewPort->getWidth(), _viewPort->getHeight());
    
    //Default transfer function
    _transferFuncTexture.init(256);
    magnet::color::TransferFunction tf;
    tf.addKnot(0,        0.91, 0.7, 0.61, 0.0);
    tf.addKnot(40.0/255, 0.91, 0.7, 0.61, 0.0);
    tf.addKnot(60.0/255, 0.91, 0.7, 0.61, 0.2);
    tf.addKnot(63.0/255, 0.91, 0.7, 0.61, 0.05);
    tf.addKnot(80.0/255, 0.91, 0.7, 0.61, 0.0);
    tf.addKnot(82.0/255, 1.0,  1.0, 0.85, 0.9);
    tf.addKnot(1.0,      1.0,  1.0, 0.85, 1.0);
    _transferFuncTexture.subImage(tf.getColorMap(), GL_RGBA);
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
	    voldata[4 * coord + 3] 
	      = inbuffer[coordCalc(x, y, z, width, height, depth)];
	  }

    _data.init(width, height, depth);
    _data.subImage(voldata, GL_RGBA);
  }
  
  void 
  RVolume::resize(size_t width, size_t height)
  {
    if (_fbo.get() != NULL) 
      _fbo->resize(width, height);
  }

  void 
  RVolume::initOpenCL() 
  {
    {
      float vertices[] = {-1,-1,-1,  1,-1,-1,  1, 1,-1, -1, 1,-1,
			  -1,-1, 1, -1, 1, 1,  1, 1, 1,  1,-1, 1};
      
      std::vector<float> VertexPos(vertices, vertices 
				   + sizeof(vertices) / sizeof(float));
      setGLPositions(VertexPos);
    }
    
    {
      int elements[] = {3,2,1,0, 6,7,1,2, 5,4,7,6, 3,0,4,5, 6,2,3,5, 7,4,0,1 };
      std::vector<int> ElementData(elements, elements + sizeof(elements) / sizeof(int));
      setGLElements(ElementData);
    }

    
  }  
  
  void 
  RVolume::glRender(magnet::GL::FBO& fbo)
  {
    if (!_visible || !_data.isValid()) return;

    //Before we render, we need the current depth buffer for depth testing.
    fbo.detach();   
    fbo.copyto(*_fbo, GL_DEPTH_BUFFER_BIT);
    fbo.attach();

    //Now bind this copied depth texture to texture unit 0
    glActiveTextureARB(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _fbo->getDepthTexture());

    _data.bind(1);
    _transferFuncTexture.bind(2);

    //Now we can render
    GLhandleARB oldshader = glGetHandleARB(GL_PROGRAM_OBJECT_ARB);

    GLfloat FocalLength = 1.0f / std::tan(_viewPort->getFOVY() * (M_PI / 360.0f));

    _shader.attach(FocalLength, _viewPort->getWidth(), 
		   _viewPort->getHeight(), _viewPort->getEyeLocation(),
		   0, 1, 2, _viewPort->getZNear(), _viewPort->getZFar(),
		   _stepSizeVal, _diffusiveLighting->get_active(),
		   _specularLighting->get_active());
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);
    RQuads::glRender();
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
    _specularLighting.reset(new Gtk::CheckButton);
    _transferFunction.reset(new magnet::gtk::TransferFunction
			    (magnet::function::MakeDelegate
			     (this, &RVolume::transferFunctionUpdated)));
    _transferFunction->set_size_request(-1, 100);
    _optList->add(*_transferFunction); _transferFunction->show();

    {//Volume renderer step size
      _stepSize.reset(new Gtk::Entry);
      Gtk::HBox* box = manage(new Gtk::HBox);	
      Gtk::Label* label = manage(new Gtk::Label("Raytrace Step Size"));
    
      box->pack_start(*label, false, false); label->show();
      box->pack_end(*_stepSize, false, false);
      _stepSize->show(); _stepSize->set_text("0.01");
      
      box->show();
      _optList->add(*box);
    }

    {//Diffusive lighting
      _diffusiveLighting.reset(new Gtk::CheckButton("Diffusive Lighting"));
      _diffusiveLighting->show();      
      _diffusiveLighting->set_active();
      _optList->add(*_diffusiveLighting);
    }

    {//Specular lighting
      _specularLighting.reset(new Gtk::CheckButton("Specular Lighting"));
      _specularLighting->show();      
      _specularLighting->set_active();
      _optList->add(*_specularLighting);
    }
    
    _optList->show();
    //Callbacks
    _stepSize->signal_changed()
      .connect(sigc::bind<Gtk::Entry&>(&magnet::Gtk::forceNumericEntry, *_stepSize));
    _stepSize->signal_activate()
      .connect(sigc::mem_fun(*this, &RVolume::guiUpdate));
    _diffusiveLighting->signal_toggled()
      .connect(sigc::mem_fun(*this, &RVolume::guiUpdate));
    _specularLighting->signal_toggled()
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

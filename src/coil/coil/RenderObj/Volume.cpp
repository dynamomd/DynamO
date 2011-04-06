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
#include <iostream>
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <boost/lexical_cast.hpp>
#include <magnet/GL/textureLoaderRAW.hpp>
#include <magnet/color/transferFunction.hpp>
#include <magnet/PNG.hpp>

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
    _data.init(width, height, depth);
    try {
      magnet::GL::loadVolumeFromRawFile(filename, _data);
    } catch (magnet::exception& tmp)
      {
	_data.deinit();
	throw;
      }
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
  RVolume::initGTK()
  {
    _optList.reset(new Gtk::VBox);//The Vbox of options   
    _specularLighting.reset(new Gtk::CheckButton);

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

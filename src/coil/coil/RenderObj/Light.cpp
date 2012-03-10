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

#include <coil/RenderObj/Light.hpp>
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <magnet/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

extern const guint8 Light_Icon[];
extern const size_t Light_Icon_size;

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RLight::getIcon()
  {
    return Gdk::Pixbuf::create_from_inline(Light_Icon_size, Light_Icon);
  }

  void 
  RLight::deinit()
  {
    _sphereShader.deinit();
    _glposition.deinit();
    _context.reset();
  }
  
  void 
  RLight::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue) 
  {
    RenderObj::init(systemQueue);

    _sphereShader.build();

    magnet::math::Vector loc = getEyeLocationObjSpace();
    GLfloat pos[3] = {loc[0], loc[1], loc[2]};
    std::vector<GLfloat> position(pos, pos + 3);
    _glposition.init(position);

    _context = magnet::GL::Context::getContext();
    initGTK();
  }

  void 
  RLight::glRender(const magnet::GL::Camera& cam, 
		   RenderMode mode)
  {
    if (!_visible) return;
    
    using namespace magnet::GL;

    if (mode & RenderObj::COLOR)
      {
	magnet::math::Vector loc = getEyeLocationObjSpace();
	GLfloat pos[3] = {loc[0], loc[1], loc[2]};
	std::vector<GLfloat> position(pos, pos + 3);
	_glposition.init(position);

	_context->cleanupAttributeArrays();
	//Set the normals to zero so it is fully illuminated
	_context->setAttribute(Context::instanceScaleAttrIndex, 0.05, 0.05, 0.05, 1);
	_context->setAttribute(Context::vertexColorAttrIndex, 1, 1, 1, 1);
	
	if (_context->testExtension("GL_ARB_sample_shading"))
	  {
#ifndef GL_SAMPLE_SHADING
# define GL_SAMPLE_SHADING GL_SAMPLE_SHADING_ARB
#endif
	    glEnable(GL_SAMPLE_SHADING);
	    glMinSampleShadingARB(1.0);
	  }

	_sphereShader.attach();
	_sphereShader["ProjectionMatrix"] = cam.getProjectionMatrix();
	_sphereShader["ViewMatrix"] = cam.getViewMatrix();
	_sphereShader["global_scale"] = GLfloat(1.0);
	_glposition.drawArray(magnet::GL::element_type::POINTS, 3);
	_sphereShader.detach();

	if (_context->testExtension("GL_ARB_sample_shading"))
	  glDisable(GL_SAMPLE_SHADING);
      }
  }

  void
  RLight::initGTK()
  {
    _optList.reset(new Gtk::VBox);

    { //Intensity
      Gtk::HBox* box = manage(new Gtk::HBox);
      box->show();
      
      {
	Gtk::Label* label = manage(new Gtk::Label("Intensity and Color", 0.95, 0.5));
	box->pack_start(*label, true, true); 
	label->show();
      }

      _intensityEntry.reset(new Gtk::Entry);
      box->pack_start(*_intensityEntry, false, false);
      _intensityEntry->show(); _intensityEntry->set_width_chars(7);
      _intensityEntry->set_text(boost::lexical_cast<std::string>(_intensity));

      _intensityEntry->signal_changed()
	.connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_intensityEntry));
      _intensityEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      {
	_lightColor.reset(new Gtk::ColorButton);
	_lightColor->set_use_alpha(false);
	Gdk::Color color;
	color.set_rgb(_color[0] * G_MAXUSHORT, _color[1] * G_MAXUSHORT, _color[2] * G_MAXUSHORT);
	_lightColor->set_color(color);
	box->pack_start(*_lightColor, false, false);
	_lightColor->show();

	_lightColor->signal_color_set().connect(sigc::mem_fun(*this, &RLight::guiUpdate));
      }

      {
	Gtk::Label* label = manage(new Gtk::Label("Attenuation", 0.95, 0.5));
	box->pack_start(*label, true, true);
	label->show();
      }

      _attenuationEntry.reset(new Gtk::Entry);
      box->pack_start(*_attenuationEntry, false, false);
      _attenuationEntry->show(); _attenuationEntry->set_width_chars(7);
      _attenuationEntry->set_text(boost::lexical_cast<std::string>(_attenuation));

      _attenuationEntry->signal_changed()
	.connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_attenuationEntry));
      _attenuationEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      _optList->pack_start(*box, false, false); 
    }

    { //Specular
      Gtk::HBox* box = manage(new Gtk::HBox);
      box->show();

      {
	Gtk::Label* label = manage(new Gtk::Label("Specular Exponent", 0.95, 0.5));
	box->pack_start(*label, true, true); 
	label->show();
      }

      _specularExponentEntry.reset(new Gtk::Entry);
      box->pack_start(*_specularExponentEntry, false, false);
      _specularExponentEntry->show(); _specularExponentEntry->set_width_chars(7);
      _specularExponentEntry->set_text(boost::lexical_cast<std::string>(_specularExponent));

      _specularExponentEntry->signal_changed()
	.connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_specularExponentEntry));
      _specularExponentEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      {
	Gtk::Label* label = manage(new Gtk::Label("Specular Strength", 0.95, 0.5));
	box->pack_start(*label, true, true); 
	label->show();
      }

      _specularFactorEntry.reset(new Gtk::Entry);
      box->pack_start(*_specularFactorEntry, false, false);
      _specularFactorEntry->show(); _specularFactorEntry->set_width_chars(7);
      _specularFactorEntry->set_text(boost::lexical_cast<std::string>(_specularFactor));

      _specularFactorEntry->signal_changed()
	.connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_specularFactorEntry));
      _specularFactorEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      _optList->pack_start(*box, false, false);
    }

    { //Specular
      Gtk::HBox* box = manage(new Gtk::HBox);
      box->show();

      {
	Gtk::Label* label = manage(new Gtk::Label("Position", 0.95, 0.5));
	box->pack_start(*label, true, true); 
	label->show();
      }

      _positionXEntry.reset(new Gtk::Entry);
      box->pack_start(*_positionXEntry, false, false);
      _positionXEntry->show();
      _positionXEntry->set_width_chars(7);
      _positionXEntry->set_text(boost::lexical_cast<std::string>(getEyeLocationObjSpace()[0]));
      _positionXEntry->signal_changed().connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_specularFactorEntry));
      _positionXEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      _positionYEntry.reset(new Gtk::Entry);
      box->pack_start(*_positionYEntry, false, false);
      _positionYEntry->show();
      _positionYEntry->set_width_chars(7);
      _positionYEntry->set_text(boost::lexical_cast<std::string>(getEyeLocationObjSpace()[1]));
      _positionYEntry->signal_changed().connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_specularFactorEntry));
      _positionYEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      _positionZEntry.reset(new Gtk::Entry);
      box->pack_start(*_positionZEntry, false, false);
      _positionZEntry->show();
      _positionZEntry->set_width_chars(7);
      _positionZEntry->set_text(boost::lexical_cast<std::string>(getEyeLocationObjSpace()[2]));
      _positionZEntry->signal_changed().connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, *_specularFactorEntry));
      _positionZEntry->signal_activate().connect(sigc::mem_fun(*this, &RLight::guiUpdate));

      _optList->pack_start(*box, false, false);
    }

    _optList->show();

    guiUpdate();
  }

  void
  RLight::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _optList->unparent();
    win->add(*_optList);
    win->show();
  }

  void 
  RLight::pickingRender(const magnet::GL::Camera& cam, 
			const uint32_t offset)
  {
    if (!_visible) return;
    
    using namespace magnet::GL;

    magnet::math::Vector loc = getEyeLocationObjSpace();
    GLfloat pos[3] = {loc[0], loc[1], loc[2]};
    std::vector<GLfloat> position(pos, pos + 3);
    _glposition.init(position);
    
    _context->cleanupAttributeArrays();
    //Set the normals to zero so it is fully illuminated
    _context->setAttribute(Context::instanceScaleAttrIndex, 0.05, 0.05, 0.05, 1);
    
    _context->setAttribute(Context::vertexColorAttrIndex, 
			   (offset % 256) / 255.0, 
			   ((offset / 256) % 256) / 255.0, 
			   (((offset / 256) / 256) % 256) / 255.0, 
			   ((((offset / 256) / 256) / 256) % 256) / 255.0);
    
    _sphereShader.attach();
    _sphereShader["ProjectionMatrix"] = cam.getProjectionMatrix();
    _sphereShader["ViewMatrix"] = cam.getViewMatrix();
    _sphereShader["global_scale"] = GLfloat(1.0);
    _glposition.drawArray(magnet::GL::element_type::POINTS, 3);
    _sphereShader.detach();
  }

  void 
  RLight::guiUpdate()
  {
    try { _intensity = boost::lexical_cast<float>(_intensityEntry->get_text()); } catch (...) {}
    try { _attenuation = boost::lexical_cast<float>(_attenuationEntry->get_text()); } catch (...) {}
    try { _specularExponent = boost::lexical_cast<float>(_specularExponentEntry->get_text()); } catch (...) {}
    try { _specularFactor = boost::lexical_cast<float>(_specularFactorEntry->get_text()); } catch (...) {}

    Gdk::Color color = _lightColor->get_color();
    _color[0] = GLfloat(color.get_red()) / G_MAXUSHORT;
    _color[1] = GLfloat(color.get_green()) / G_MAXUSHORT;
    _color[2] = GLfloat(color.get_blue()) / G_MAXUSHORT;

    try {
      magnet::math::Vector vec;
      vec[0] = boost::lexical_cast<float>(_positionXEntry->get_text());
      vec[1] = boost::lexical_cast<float>(_positionYEntry->get_text());
      vec[2] = boost::lexical_cast<float>(_positionZEntry->get_text());
      setEyeLocationObjSpace(vec);
    } catch (...) {}
  }
}

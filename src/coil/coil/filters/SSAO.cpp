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

#include <coil/filters/SSAO.hpp>

extern const char _binary_src_coil_coil_filters_SSAOwin_gladexml_start[];
extern const char _binary_src_coil_coil_filters_SSAOwin_gladexml_end[];

namespace coil 
{
  SSAOWrapper::SSAOWrapper()
  { 
    _radius = 0.03;
    _totStrength = 10000.0;
    _strength = 0.25;
    _offset = 13.11;

    _filter.build(); 

    Glib::ustring glade_data
      (reinterpret_cast<const char *>(_binary_src_coil_coil_filters_SSAOwin_gladexml_start), 
       _binary_src_coil_coil_filters_SSAOwin_gladexml_end
       -_binary_src_coil_coil_filters_SSAOwin_gladexml_start);
    
    _refXml = Gtk::Builder::create_from_string(glade_data);

    {
      Gtk::SpinButton* btn;

      _refXml->get_widget("radiusSetting", btn);
      btn->set_value(_radius);
      btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));

      _refXml->get_widget("totStrengthSetting", btn);
      btn->set_value(_totStrength);
      btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));

      _refXml->get_widget("strengthSetting", btn);
      btn->set_value(_strength);
      btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));

      _refXml->get_widget("offsetSetting", btn);
      btn->set_value(_offset);
      btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    }

    glGenTextures( 1, &_randomTexture);
    glBindTexture(GL_TEXTURE_2D, _randomTexture);

    const size_t randomTextureSize = 64;

    std::vector<GLubyte> texture;
    texture.resize(3 * randomTextureSize * randomTextureSize);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    srand(120121);
    for (size_t i(0); i < 3 * randomTextureSize * randomTextureSize; ++i)
      texture[i] = (rand() * 255.0) / RAND_MAX;
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, randomTextureSize, randomTextureSize, 
		 0, GL_RGB, GL_UNSIGNED_BYTE, &texture[0]);
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
  }

  SSAOWrapper::~SSAOWrapper()
  { 
    Gtk::Window* win;
    _refXml->get_widget("SSAOeditWindow", win);
    win->hide();

    glDeleteTextures(1, &_randomTexture);
  }

  void SSAOWrapper::edit()
  {
    Gtk::Window* win;
    _refXml->get_widget("SSAOeditWindow", win);
    win->show();
  }

  void SSAOWrapper::settingsCallback()
  {
    Gtk::SpinButton* btn;
    _refXml->get_widget("radiusSetting", btn);
    _radius = btn->get_value();

    _refXml->get_widget("totStrengthSetting", btn);
    _totStrength = btn->get_value();

    _refXml->get_widget("strengthSetting", btn);
    _strength = btn->get_value();

    _refXml->get_widget("offsetSetting", btn);
    _offset = btn->get_value();
  }

  void SSAOWrapper::invoke(GLuint colorTextureUnit, 
			   size_t width, size_t height) 
  {
    glActiveTextureARB(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_2D, _randomTexture);

    _filter.invoke(colorTextureUnit, 1, 2, 7, width, height, 
		   _radius, _totStrength, _strength, _offset); 
  }

}


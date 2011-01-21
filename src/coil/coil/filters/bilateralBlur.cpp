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

#include <coil/filters/bilateralBlur.hpp>

extern const char _binary_src_coil_coil_filters_bilateralBlurWin_gladexml_start[];
extern const char _binary_src_coil_coil_filters_bilateralBlurWin_gladexml_end[];

namespace coil 
{
  ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////Blur Phase/////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  BilateralBlurWrapper::BilateralBlurWrapper()
  {
    _radius = 1;
    _totStrength = 0.01245;

    _filter.build(); 

    Glib::ustring glade_data
      (reinterpret_cast<const char *>(_binary_src_coil_coil_filters_bilateralBlurWin_gladexml_start), 
       _binary_src_coil_coil_filters_bilateralBlurWin_gladexml_end
       -_binary_src_coil_coil_filters_bilateralBlurWin_gladexml_start);
    
    _refXml = Gtk::Builder::create_from_string(glade_data);

    {
      Gtk::HScale* btn;
      _refXml->get_widget("radiusSetting", btn);
      btn->set_value(_radius);
      btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));

      _refXml->get_widget("zdiffsetting", btn);
      btn->set_range(0.0001, 1.0);
      btn->set_value(_totStrength); 
     btn->signal_value_changed()
	.connect(sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));
    }

    glGenTextures( 1, &_randomTexture);
    glBindTexture(GL_TEXTURE_2D, _randomTexture);

    std::vector<GLubyte> texture;
    texture.resize(3 * 64 * 64);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    srand(120121);
    for (size_t i(0); i < 3 * 64 * 64; ++i)
      texture[i] = (rand() * 255.0) / RAND_MAX;
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 64, 64, 0, GL_RGB, GL_UNSIGNED_BYTE, &texture[0]);
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
  }

  BilateralBlurWrapper::~BilateralBlurWrapper()
  { 
    Gtk::Window* win;
    _refXml->get_widget("BilateralBlurWindow", win);
    win->hide();

    glDeleteTextures(1, &_randomTexture);
  }

  void BilateralBlurWrapper::edit()
  {
    Gtk::Window* win;
    _refXml->get_widget("BilateralBlurWindow", win);
    win->show();
  }

  void BilateralBlurWrapper::settingsCallback()
  {
    Gtk::HScale* btn;
    _refXml->get_widget("radiusSetting", btn);
    _radius = btn->get_value();

    _refXml->get_widget("zdiffsetting", btn);
    _totStrength = btn->get_value();
  }

  void BilateralBlurWrapper::invoke(GLuint colorTextureUnit, 
				    size_t width, size_t height) 
  {
    _filter.invoke(colorTextureUnit, 2, width, height, _radius, _totStrength); 
  }
}

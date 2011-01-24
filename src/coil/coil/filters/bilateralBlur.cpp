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

    //Build the controls

    {
      Gtk::VBox* labelVbox = manage(new Gtk::VBox);
      Gtk::Label* Label1 = manage(new Gtk::Label("Radius"));
      Gtk::Label* Label2 = manage(new Gtk::Label("Depth Cutoff"));
      labelVbox->add(*Label1); Label1->show();
      labelVbox->add(*Label2); Label2->show();
      _optlist.add(*labelVbox);
      labelVbox->show();
    }

    {
      Gtk::VBox* sliderVbox = manage(new Gtk::VBox);
      sliderVbox->add(_radiusSlider);
      sliderVbox->add(_zdiffSlider);
      _optlist.add(*sliderVbox);
      sliderVbox->show();
    }
    _optlist.show();

    _radiusSlider.set_range(1,20);
    _radiusSlider.set_increments(1,1);
    _radiusSlider.set_value(_radius);
    _radiusSlider.signal_value_changed()
      .connect(sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));
    _radiusSlider.show();

    _zdiffSlider.set_increments(0.0001,0.0001);
    _zdiffSlider.set_range(0.0001, 1.0);
    _zdiffSlider.set_digits(4);
    _zdiffSlider.set_value(_totStrength); 
    _zdiffSlider.signal_value_changed()
      .connect(sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));
    _zdiffSlider.show();
  }

  BilateralBlurWrapper::~BilateralBlurWrapper()
  { 
    glDeleteTextures(1, &_randomTexture);
  }

  void BilateralBlurWrapper::showControls(Gtk::ScrolledWindow* start)
  {
    _optlist.unparent();
    start->add(_optlist);
    start->show();


  }

  void BilateralBlurWrapper::settingsCallback()
  {
    _radius = _radiusSlider.get_value();
    _totStrength = _zdiffSlider.get_value();
  }

  void BilateralBlurWrapper::invoke(GLuint colorTextureUnit, 
				    size_t width, size_t height) 
  {
    _filter.invoke(colorTextureUnit, 2, width, height, _radius, _totStrength); 
  }
}

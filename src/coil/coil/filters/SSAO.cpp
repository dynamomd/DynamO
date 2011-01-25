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

namespace coil 
{
  SSAOWrapper::SSAOWrapper()
  { 
    _radius = 0.03;
    _totStrength = 10;
    _strength = 0.25;
    _offset = 13.11;

    _filter.build(); 

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

    //Build the controls
    {
      Gtk::VBox* labelVbox = manage(new Gtk::VBox);
      Gtk::Label* Label1 = manage(new Gtk::Label("Radius"));
      Gtk::Label* Label2 = manage(new Gtk::Label("Magnitude"));
      Gtk::Label* Label3 = manage(new Gtk::Label("Drop off"));
      Gtk::Label* Label4 = manage(new Gtk::Label("Random number offset"));
      labelVbox->add(*Label1); Label1->show();
      labelVbox->add(*Label2); Label2->show();
      labelVbox->add(*Label3); Label3->show();
      labelVbox->add(*Label4); Label4->show();
      _optlist.add(*labelVbox);
      labelVbox->show();
    }

    {
      Gtk::VBox* sliderVbox = manage(new Gtk::VBox);
      sliderVbox->add(_radiusSlider);
      sliderVbox->add(_totStrengthSlider);
      sliderVbox->add(_strengthSlider);
      sliderVbox->add(_randomOffsetSlider);

      _optlist.add(*sliderVbox);
      sliderVbox->show();
    }
    _optlist.show();

    _radiusSlider.set_range(0.001,1);
    _radiusSlider.set_increments(1,1);
    _radiusSlider.set_digits(3);
    _radiusSlider.set_value(_radius);
    _radiusSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _radiusSlider.show();

    _totStrengthSlider.set_range(0.1,10);
    _totStrengthSlider.set_increments(1,1);
    _totStrengthSlider.set_digits(1);
    _totStrengthSlider.set_value(_totStrength);
    _totStrengthSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _totStrengthSlider.show();

    _strengthSlider.set_range(0.01,1);
    _strengthSlider.set_increments(0.1,0.1);
    _strengthSlider.set_digits(2);
    _strengthSlider.set_value(_strength);
    _strengthSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _strengthSlider.show();

    _randomOffsetSlider.set_range(1,100);
    _randomOffsetSlider.set_increments(0.1,0.1);
    _randomOffsetSlider.set_digits(2);
    _randomOffsetSlider.set_value(_offset);
    _randomOffsetSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _randomOffsetSlider.show();
  }

  void 
  SSAOWrapper::showControls(Gtk::ScrolledWindow* start)
  {
    _optlist.unparent();
    start->add(_optlist);
    start->show();
  }
  
  void 
  SSAOWrapper::settingsCallback()
  {
    _radius = _radiusSlider.get_value();
    _totStrength = _totStrengthSlider.get_value();
    _strength = _strengthSlider.get_value();
    _offset = _randomOffsetSlider.get_value();
  }


  SSAOWrapper::~SSAOWrapper()
  { 
    glDeleteTextures(1, &_randomTexture);
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


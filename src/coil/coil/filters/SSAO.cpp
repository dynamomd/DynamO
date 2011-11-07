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

#include <gtkmm.h>
#include <coil/filters/SSAO.hpp>

namespace coil 
{
  SSAOWrapper::SSAOWrapper()
  { 
    _radius = 0.005;
    _totStrength = 1;
    _dropoff = 0.05;

    _filter.build(); 

    glGenTextures( 1, &_randomTexture);
    glBindTexture(GL_TEXTURE_2D, _randomTexture);

    std::vector<GLubyte> texture;
    texture.resize(3 * _randomTextureSize * _randomTextureSize);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    srand(120121);
    for (size_t i(0); i < 3 * _randomTextureSize * _randomTextureSize; ++i)
      texture[i] = (rand() * 255.0) / RAND_MAX;
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, _randomTextureSize, _randomTextureSize, 
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
      labelVbox->add(*Label1); Label1->show();
      labelVbox->add(*Label2); Label2->show();
      labelVbox->add(*Label3); Label3->show();
      _optlist.add(*labelVbox);
      labelVbox->show();
    }

    {
      Gtk::VBox* sliderVbox = manage(new Gtk::VBox);
      sliderVbox->add(_radiusSlider);
      sliderVbox->add(_totStrengthSlider);
      sliderVbox->add(_dropoffSlider);

      _optlist.add(*sliderVbox);
      sliderVbox->show();
    }
    _optlist.show();

    _radiusSlider.set_range(0.00001, 0.01);
    _radiusSlider.set_increments(1,1);
    _radiusSlider.set_digits(6);
    _radiusSlider.set_value(_radius);
    _radiusSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _radiusSlider.show();

    _totStrengthSlider.set_range(0,2);
    _totStrengthSlider.set_increments(1,1);
    _totStrengthSlider.set_digits(3);
    _totStrengthSlider.set_value(_totStrength);
    _totStrengthSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _totStrengthSlider.show();

    _dropoffSlider.set_range(0.00001,0.20);
    _dropoffSlider.set_increments(0.1,0.1);
    _dropoffSlider.set_digits(5);
    _dropoffSlider.set_value(_dropoff);
    _dropoffSlider.signal_value_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _dropoffSlider.show();
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
    _dropoff = _dropoffSlider.get_value();
  }


  SSAOWrapper::~SSAOWrapper()
  { 
    glDeleteTextures(1, &_randomTexture);
  }

  void SSAOWrapper::invoke(GLint colorTextureUnit, 
			   size_t width, size_t height,
			   const magnet::GL::Camera& vp)
  {
    glActiveTextureARB(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_2D, _randomTexture);

    _filter.attach();
    _filter["radius"] = _radius;
    _filter["totStrength"] = _totStrength;
    _filter["depthDropoff"] = _dropoff;
    _filter["offset"] = GLfloat(std::max(width, height)) / _randomTextureSize;
    _filter["NormalsTex"] = 1;
    _filter["EyePosTex"] = 2;
    _filter["rnm"] = 7;
    _filter.invoke();
    _filter.detach();
  }
}


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
#include <magnet/gtk/numericEntry.hpp>
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

    _radiusSlider.set_text(boost::lexical_cast<std::string>(_radius));
    _radiusSlider.signal_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _radiusSlider.show();

    _totStrengthSlider.set_text(boost::lexical_cast<std::string>(_totStrength));
    _totStrengthSlider.signal_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
    _totStrengthSlider.show();
    
    _dropoffSlider.set_text(boost::lexical_cast<std::string>(_dropoff));
    _dropoffSlider.signal_changed().connect(sigc::mem_fun(this, &SSAOWrapper::settingsCallback));
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
    magnet::gtk::forceNumericEntry(_radiusSlider);
    magnet::gtk::forceNumericEntry(_totStrengthSlider);
    magnet::gtk::forceNumericEntry(_dropoffSlider);

    try { _radius = boost::lexical_cast<float>(_radiusSlider.get_text()); } catch (...) {};
    try { _totStrength = boost::lexical_cast<float>(_totStrengthSlider.get_text()); } catch (...) {};
    try { _dropoff = boost::lexical_cast<float>(_dropoffSlider.get_text()); } catch (...) {};
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


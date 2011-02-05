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

#pragma once
#include "filter.hpp"
#include <magnet/GL/DOF.hpp>

namespace coil 
{
  class DOFFilter: public filter
  {
  public:
    DOFFilter():
      _focalLength(0.025f), _focalWidth(0.01f)
    { 
      _filter.build();
      
      //Build the controls
      {
	Gtk::VBox* labelVbox = manage(new Gtk::VBox);
	Gtk::Label* Label1 = manage(new Gtk::Label("Focal Length"));
	Gtk::Label* Label2 = manage(new Gtk::Label("Focal Width"));
	labelVbox->add(*Label1); Label1->show();
	labelVbox->add(*Label2); Label2->show();
	_optlist.add(*labelVbox);
	labelVbox->show();
      }
      
      {
	Gtk::VBox* sliderVbox = manage(new Gtk::VBox);
	sliderVbox->add(_focalLengthSlider);
	sliderVbox->add(_focalWidthSlider);
	_optlist.add(*sliderVbox);
	sliderVbox->show();
      }
      _optlist.show();
      
      _focalLengthSlider.set_range(0, 0.1);
      _focalLengthSlider.set_increments(0.01,0.01);
      _focalLengthSlider.set_digits(3);
      _focalLengthSlider.set_value(_focalLength);
      _focalLengthSlider.signal_value_changed().connect(sigc::mem_fun(this, &DOFFilter::settingsCallback));
      _focalLengthSlider.show();
      
      _focalWidthSlider.set_increments(0.001,0.001);
      _focalWidthSlider.set_range(0.0001, 0.25);
      _focalWidthSlider.set_digits(4);
      _focalWidthSlider.set_value(_focalWidth); 
      _focalWidthSlider.signal_value_changed().connect(sigc::mem_fun(this, &DOFFilter::settingsCallback));
      _focalWidthSlider.show();
    }

    inline virtual size_t type_id() { return detail::filterEnum<DOFFilter>::val; }
    inline virtual void invoke(GLuint colorTextureUnit, size_t width, size_t height,
			       const magnet::GL::viewPort& vp)
    { _filter.invoke(colorTextureUnit, 0, 2, _focalLength, _focalWidth, width, height); }

    inline virtual bool needsNormalDepth()  { return false; }

    virtual void showControls(Gtk::ScrolledWindow* start)
    {
      _optlist.unparent();
      start->add(_optlist);
      start->show();
    }
  protected:
    magnet::GL::DOF _filter;

    Gtk::HScale _focalLengthSlider;
    Gtk::HScale _focalWidthSlider;
    Gtk::HBox _optlist;

    GLfloat _focalLength, _focalWidth;

    void settingsCallback()
    {
      _focalLength = _focalLengthSlider.get_value();
      _focalWidth = _focalWidthSlider.get_value();
    }
  };
}

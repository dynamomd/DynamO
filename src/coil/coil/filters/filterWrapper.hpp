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

#pragma once

#include <gtkmm.h>
#include "filter.hpp"

namespace coil 
{
  template <class T, bool reqNormalDepth>
  class magnetFilterWrapper: public filter
  {
  public:
    magnetFilterWrapper():
      _radius(1)
    { 
      _filter.build(); 
      
      {
	Gtk::Label* Label1 = manage(new Gtk::Label("Radius"));
	_optlist.add(*Label1); Label1->show();
      }

      _radiusSlider.set_range(1, 20);
      _radiusSlider.set_increments(1,1);
      _radiusSlider.set_digits(0);
      _radiusSlider.set_value(_radius);
      _radiusSlider.signal_value_changed()
	.connect(sigc::mem_fun(this, &magnetFilterWrapper<T,reqNormalDepth>::settingsCallback));
      _optlist.add(_radiusSlider); _radiusSlider.show();
      
      _optlist.show();
    }

    inline virtual size_t type_id() { return detail::filterEnum<magnetFilterWrapper<T,reqNormalDepth> >::val; }    


    inline virtual void invoke(GLint colorTextureUnit, size_t width, size_t height, const magnet::GL::Camera&)
    { 
      _filter["u_Texture0"] = colorTextureUnit;

      std::tr1::array<GLfloat, 2> arg = {{GLfloat(_radius) / width, GLfloat(_radius) / height}};
      _filter["u_Scale"] = arg;

      _filter.invoke(); 
    }

    inline virtual bool needsNormalDepth() { return reqNormalDepth; }

    virtual void showControls(Gtk::ScrolledWindow* start)
    {
      _optlist.unparent();
      start->add(_optlist);
      start->show();
    }
    
  protected:
    T _filter;

    GLuint _radius;

    Gtk::HScale _radiusSlider;
    Gtk::HBox _optlist;

    void settingsCallback()
    {
      _radius = _radiusSlider.get_value();
    }
  };
}

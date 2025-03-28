/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
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
#include <coil/filters/filter.hpp>
#include <magnet/GL/shader/DOF.hpp>
#include <magnet/gtk/numericEntry.hpp>

namespace coil {
class DOFFilter : public Filter {
public:
  DOFFilter() : _focalLength(0.0f), _focalWidth(1.5f) {
    _filter.build();

    // Build the controls
    {
      Gtk::VBox *labelVbox = manage(new Gtk::VBox);
      Gtk::Label *Label1 = manage(new Gtk::Label("Focal Length (0=auto)"));
      Gtk::Label *Label2 = manage(new Gtk::Label("Focal Width"));
      labelVbox->add(*Label1);
      Label1->show();
      labelVbox->add(*Label2);
      Label2->show();
      _optlist.add(*labelVbox);
      labelVbox->show();
    }

    {
      Gtk::VBox *sliderVbox = manage(new Gtk::VBox);
      sliderVbox->add(_focalLengthSlider);
      sliderVbox->add(_focalWidthSlider);
      _optlist.add(*sliderVbox);
      sliderVbox->show();
    }
    _optlist.show();

    _focalLengthSlider.set_text("0.0");
    _focalWidthSlider.set_text(boost::lexical_cast<std::string>(_focalWidth));

    _focalLengthSlider.show();
    _focalWidthSlider.show();

    _focalLengthSlider.signal_changed().connect(
        sigc::bind(&magnet::gtk::forceNumericEntry, &_focalLengthSlider));

    _focalWidthSlider.signal_changed().connect(
        sigc::bind(&magnet::gtk::forceNumericEntry, &_focalWidthSlider));

    _focalLengthSlider.signal_activate().connect(
        sigc::mem_fun(this, &DOFFilter::settingsCallback));
    _focalWidthSlider.signal_activate().connect(
        sigc::mem_fun(this, &DOFFilter::settingsCallback));
  }

  inline virtual size_t type_id() { return detail::filterEnum<DOFFilter>::val; }

  inline virtual void invoke(GLint colorTextureUnit, size_t width,
                             size_t height, const magnet::GL::Camera &vp) {
    _filter.attach();
    _filter["u_Texture0"] = colorTextureUnit;
    _filter["u_Texture1"] = 0;
    _filter["u_Texture2"] = 2;
    _filter["focalDistance"] = _focalLength;
    _filter["focalRange"] = _focalWidth;
    _filter.invoke();
    _filter.detach();
  }

  virtual void showControls(Gtk::ScrolledWindow *start) {
    _optlist.unparent();
    start->add(_optlist);
    start->show();
  }

protected:
  magnet::GL::shader::DOFShader _filter;

  Gtk::Entry _focalLengthSlider;
  Gtk::Entry _focalWidthSlider;
  Gtk::HBox _optlist;

  GLfloat _focalLength, _focalWidth;

  void settingsCallback() {
    try {
      _focalLength = boost::lexical_cast<double>(_focalLengthSlider.get_text());
    } catch (...) {
    }
    try {
      _focalWidth = boost::lexical_cast<double>(_focalWidthSlider.get_text());
    } catch (...) {
    }
  }
};
} // namespace coil

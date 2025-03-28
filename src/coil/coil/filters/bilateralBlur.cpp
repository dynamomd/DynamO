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

#include <coil/filters/bilateralBlur.hpp>
#include <magnet/gtk/numericEntry.hpp>

namespace coil {
///////////////////////////////////////////////////////////////////////////
////////////////////////////Blur Phase/////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

BilateralBlurWrapper::BilateralBlurWrapper() : _radius(1), _zdiff(0.01) {
  _filter.build();

  // Build the controls
  {
    Gtk::VBox *labelVbox = manage(new Gtk::VBox);
    Gtk::Label *Label1 = manage(new Gtk::Label("Radius"));
    Gtk::Label *Label2 = manage(new Gtk::Label("Depth Cutoff"));
    labelVbox->add(*Label1);
    Label1->show();
    labelVbox->add(*Label2);
    Label2->show();
    _optlist.add(*labelVbox);
    labelVbox->show();
  }

  {
    Gtk::VBox *sliderVbox = manage(new Gtk::VBox);
    sliderVbox->add(_radiusSlider);
    sliderVbox->add(_zdiffEntry);
    _optlist.add(*sliderVbox);
    sliderVbox->show();
  }
  _optlist.show();

  _radiusSlider.set_range(1, 20);
  _radiusSlider.set_increments(1, 1);
  _radiusSlider.set_digits(0);
  _radiusSlider.set_value(_radius);
  _radiusSlider.signal_value_changed().connect(
      sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));
  _radiusSlider.show();

  _zdiffEntry.set_text(boost::lexical_cast<std::string>(_zdiff));

  _zdiffEntry.signal_changed().connect(
      sigc::bind(&magnet::gtk::forceNumericEntry, &_zdiffEntry));

  _zdiffEntry.signal_activate().connect(
      sigc::mem_fun(this, &BilateralBlurWrapper::settingsCallback));

  _zdiffEntry.show();
}

void BilateralBlurWrapper::showControls(Gtk::ScrolledWindow *start) {
  _optlist.unparent();
  start->add(_optlist);
  start->show();
}

void BilateralBlurWrapper::settingsCallback() {
  _radius = _radiusSlider.get_value();

  try {
    _zdiff = boost::lexical_cast<double>(_zdiffEntry.get_text());
  } catch (...) {
  }
}

void BilateralBlurWrapper::invoke(GLint colorTextureUnit, size_t width,
                                  size_t height, const magnet::GL::Camera &vp) {
  _filter.attach();
  _filter["totStrength"] = _zdiff;
  _filter["radius"] = _radius;
  _filter["ImageTex"] = colorTextureUnit;
  _filter["EyePosTex"] = 3;
  _filter.invoke();
  _filter.detach();
}
} // namespace coil

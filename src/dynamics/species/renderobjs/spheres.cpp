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

#ifdef DYNAMO_visualizer
# include "spheres.hpp"
# include <magnet/gtk/numericEntry.hpp>
# include <boost/lexical_cast.hpp>

SphereParticleRenderer::SphereParticleRenderer(size_t N, std::string name, 
					       magnet::function::Delegate1<magnet::CL::CLGLState&, void> updateColorFunc):
  RTSpheres(N, name),
    _mode(COLOR_BY_ID),
  _updateColorFunc(updateColorFunc)
{
  for (size_t i(0); i < 4; ++i)
    _colorStatic.s[i] = _colorFixed.s[i] = 255;
}

void 
SphereParticleRenderer::initGTK()
{
  _optList.reset(new Gtk::VBox);//The Vbox of options
  _singleColorMode.reset(new Gtk::RadioButton("Single Color"));
  _colorByIDMode.reset(new Gtk::RadioButton("Color by ID"));
  _colorBySpeedMode.reset(new Gtk::RadioButton("Color by Speed"));
    
  _RFixed.reset(new Gtk::SpinButton);
  _GFixed.reset(new Gtk::SpinButton);
  _BFixed.reset(new Gtk::SpinButton);
  _AFixed.reset(new Gtk::SpinButton);
  _characteristicV.reset(new Gtk::Entry);
    
  _colorIfStatic.reset(new Gtk::CheckButton("Color if Static/Sleeping"));
  _RStatic.reset(new Gtk::SpinButton);
  _GStatic.reset(new Gtk::SpinButton);
  _BStatic.reset(new Gtk::SpinButton);
  _AStatic.reset(new Gtk::SpinButton);
    
    
  {//Static color and RGBA boxes
    Gtk::HBox* box = manage(new Gtk::HBox);
      
    box->pack_start(*_colorIfStatic, true, true);_colorIfStatic->show();
    _colorIfStatic->set_active();
      
    Gtk::Label* label = manage(new Gtk::Label("RGBA"));
    box->pack_start(*label, false, false); label->show();
      
    _RStatic->set_increments(1.0, 1.0);
    _RStatic->set_range(0.0, 255.0);
    _RStatic->set_value(_colorStatic.s[0]);
    _GStatic->set_increments(1.0, 1.0);
    _GStatic->set_range(0.0, 255.0);
    _GStatic->set_value(_colorStatic.s[1]);
    _BStatic->set_increments(1.0, 1.0);
    _BStatic->set_range(0.0, 255.0);
    _BStatic->set_value(_colorStatic.s[2]);
    _AStatic->set_increments(1.0, 1.0);
    _AStatic->set_range(0.0, 255.0);
    _AStatic->set_value(_colorStatic.s[3]);      
      
    box->pack_start(*_RStatic, false, false); _RStatic->show();
    box->pack_start(*_GStatic, false, false); _GStatic->show();
    box->pack_start(*_BStatic, false, false); _BStatic->show();
    box->pack_start(*_AStatic, false, false); _AStatic->show();
      
    _optList->add(*box);
    box->show();
      
    //Horizontal seperator
    Gtk::HSeparator* line = manage(new Gtk::HSeparator);
    line->show();
    _optList->add(*line);
  }
    
  {//Single color and RGBA boxes
    Gtk::HBox* box = manage(new Gtk::HBox);
      
    box->pack_start(*_singleColorMode, true, true);_singleColorMode->show();
    if (_mode == SINGLE_COLOR) _singleColorMode->set_active();
      
    Gtk::Label* label = manage(new Gtk::Label("RGBA"));
    box->pack_start(*label, false, false); label->show();
      
    _RFixed->set_increments(1.0, 1.0);
    _RFixed->set_range(0.0, 255.0);
    _RFixed->set_value(_colorFixed.s[0]);
      
    _GFixed->set_increments(1.0, 1.0);
    _GFixed->set_range(0.0, 255.0);
    _GFixed->set_value(_colorFixed.s[1]);
      
    _BFixed->set_increments(1.0, 1.0);
    _BFixed->set_range(0.0, 255.0);
    _BFixed->set_value(_colorFixed.s[2]);
      
    _AFixed->set_increments(1.0, 1.0);
    _AFixed->set_range(0.0, 255.0);
    _AFixed->set_value(_colorFixed.s[3]);      
      
    box->pack_start(*_RFixed, false, false); _RFixed->show();
    box->pack_start(*_GFixed, false, false); _GFixed->show();
    box->pack_start(*_BFixed, false, false); _BFixed->show();
    box->pack_start(*_AFixed, false, false); _AFixed->show();
      
    _optList->add(*box);
    box->show();
      
    //Horizontal seperator
    Gtk::HSeparator* line = manage(new Gtk::HSeparator);
    line->show();
    _optList->add(*line);
  }
    
  Gtk::RadioButton::Group group = _singleColorMode->get_group();
    
  {//Color by ID
    _colorByIDMode->set_group(group);
    _colorByIDMode->show();
    if (_mode == COLOR_BY_ID) _colorByIDMode->set_active();
      
    _optList->add(*_colorByIDMode);
      
    //Horizontal seperator
    Gtk::HSeparator* line = manage(new Gtk::HSeparator);
    line->show();
    _optList->add(*line);
  }
    
  {//Color by Speed
    Gtk::HBox* box = manage(new Gtk::HBox);	
      
    _colorBySpeedMode->set_group(group);
    _colorBySpeedMode->show();
    if (_mode == COLOR_BY_SPEED) _colorBySpeedMode->set_active();
      
    box->pack_start(*_colorBySpeedMode, true, true);
      
    Gtk::Label* label = manage(new Gtk::Label("Max Speed"));
    box->pack_start(*label, false, false); label->show();
      
    box->pack_start(*_characteristicV, false, false);
    _characteristicV->show();
    _characteristicV->set_text("1.0");
      
    box->show();
    _optList->add(*box);
  }
    
  _singleColorMode->signal_toggled()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _RFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _GFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _BFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _AFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
    
  _colorByIDMode->signal_toggled()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));

  _colorBySpeedMode->signal_toggled()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _characteristicV->signal_changed()
    .connect(sigc::bind<Gtk::Entry&>(&magnet::Gtk::forceNumericEntry, *_characteristicV));
  _characteristicV->signal_activate()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));


  _colorIfStatic->signal_toggled()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _RStatic->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _GStatic->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _BStatic->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));
  _AStatic->signal_value_changed()
    .connect(sigc::mem_fun(*this, &SphereParticleRenderer::guiUpdate));

  _optList->show();

  guiUpdate();
}

void 
SphereParticleRenderer::showControls(Gtk::ScrolledWindow* win)
{
  win->remove();
  _optList->unparent();
  win->add(*_optList);
  win->show();
}

void 
SphereParticleRenderer::guiUpdate()
{
  std::string val = _characteristicV->get_text();
  if (val.empty()) {val = "1.0"; _characteristicV->set_text("1.0"); }

  _scaleV = boost::lexical_cast<double>(val);
  if (_scaleV == 0) {_scaleV = 1; _characteristicV->set_text("1.0"); }
      
  if (_singleColorMode->get_active())  { _mode = SINGLE_COLOR; _recolorOnUpdate = false; }
  if (_colorByIDMode->get_active())    { _mode = COLOR_BY_ID; _recolorOnUpdate = false; }
  if (_colorBySpeedMode->get_active()) { _mode = COLOR_BY_SPEED; _recolorOnUpdate = true; }
  _colorFixed.s[0] = _RFixed->get_value();
  _colorFixed.s[1] = _GFixed->get_value();
  _colorFixed.s[2] = _BFixed->get_value();
  _colorFixed.s[3] = _AFixed->get_value();

  _colorStatic.s[0] = _RStatic->get_value();
  _colorStatic.s[1] = _GStatic->get_value();
  _colorStatic.s[2] = _BStatic->get_value();
  _colorStatic.s[3] = _AStatic->get_value();
      
  _colorStaticParticles = _colorIfStatic->get_active();

  _systemQueue->queueTask(new magnet::function::Task1<void, magnet::CL::CLGLState&>
			  (_updateColorFunc, *_CLState));
}

#endif

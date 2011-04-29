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

#include <coil/RenderObj/SphericalParticles.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <boost/lexical_cast.hpp>

RSphericalParticles::RSphericalParticles(size_t N, std::string name, size_t spheresPerObject):
  RTSpheres(N, name),
  _spheresPerObject(spheresPerObject),
  _mode(SINGLE_COLOR)
{
  _particleData.resize(N);
  _particleColorData.resize(N);

  for (size_t i(0); i < 4; ++i)
    _colorFixed.s[i] = 255;

  _colorFixed.s[0] = 255;
  _colorFixed.s[1] = 0;
  _colorFixed.s[2] = 0;
  _colorFixed.s[3] = 255;
}

void 
RSphericalParticles::sendRenderData(magnet::CL::CLGLState& CLState)
{
  CLState.getCommandQueue().enqueueWriteBuffer
    (getSphereDataBuffer(), false, 0, _N * sizeof(cl_float4), &_particleData[0]);
}

void 
RSphericalParticles::sendColorData(magnet::CL::CLGLState& CLState)
{
  const size_t segmentSize = _N / _spheresPerObject;
  for (size_t i(0); i < _spheresPerObject; ++i)
    CLState.getCommandQueue().enqueueWriteBuffer
      (getColorDataBuffer(), false, i * segmentSize * sizeof(cl_uchar4), segmentSize * sizeof(cl_uchar4), &_particleColorData[0]);
}

void 
RSphericalParticles::updateColorData(magnet::CL::CLGLState& CLState)
{
  switch (getDrawMode())
    {
    case SINGLE_COLOR:
      {
	for (size_t i(0); i < _N; ++i)
	  for (size_t cc(0); cc < 4; ++cc)
	    _particleColorData[i].s[cc] = _colorFixed.s[cc];
	break;
      }
    case COLOR_BY_ID:
      {
	for (size_t i(0); i < _N; ++i)
	  map(_particleColorData[i], ((float)(i)) / _N);
	break;
      }
    default:
      M_throw() << "Not Implemented";
    }

  _coil.getInstance().getTaskQueue()
    .queueTask(magnet::function::Task::makeTask(&RSphericalParticles::sendColorData, 
						this, CLState));
}

void 
RSphericalParticles::initGTK()
{
  _optList.reset(new Gtk::VBox);//The Vbox of options
  _colorMap.reset(new magnet::Gtk::ColorMapSelector);

  _singleColorMode.reset(new Gtk::RadioButton("Single Color"));
  _colorByIDMode.reset(new Gtk::RadioButton("Color by ID"));
    
  _RFixed.reset(new Gtk::SpinButton);
  _GFixed.reset(new Gtk::SpinButton);
  _BFixed.reset(new Gtk::SpinButton);
  _AFixed.reset(new Gtk::SpinButton);
    
  {
    _optList->add(*_colorMap); _colorMap->show();
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
    
  _optList->show();
    
  _singleColorMode->signal_toggled()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
  _colorMap->signal_changed()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
  _RFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
  _GFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
  _BFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
  _AFixed->signal_value_changed()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));
    
  _colorByIDMode->signal_toggled()
    .connect(sigc::mem_fun(*this, &RSphericalParticles::guiUpdate));

  guiUpdate();
}

void 
RSphericalParticles::showControls(Gtk::ScrolledWindow* win)
{
  win->remove();
  _optList->unparent();
  win->add(*_optList);
  win->show();
}

void 
RSphericalParticles::guiUpdate()
{
  if (_singleColorMode->get_active())  { _mode = SINGLE_COLOR; _recolorOnUpdate = false; }
  if (_colorByIDMode->get_active())    { _mode = COLOR_BY_ID; _recolorOnUpdate = false; }
  _colorFixed.s[0] = _RFixed->get_value();
  _colorFixed.s[1] = _GFixed->get_value();
  _colorFixed.s[2] = _BFixed->get_value();
  _colorFixed.s[3] = _AFixed->get_value();

  //! Get the 
  _systemQueue->queueTask(magnet::function::Task::makeTask
			  (&RSphericalParticles::updateColorData, this, *_CLState));
}

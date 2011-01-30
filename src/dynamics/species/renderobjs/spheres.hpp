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
#ifdef DYNAMO_visualizer
# include <gtkmm.h>
# include <coil/RenderObj/Spheres.hpp>
# include <memory>

class SphereParticleRenderer: public RTSpheres
{
public:
  SphereParticleRenderer(size_t N, std::string name, 
			 magnet::function::Delegate1<magnet::CL::CLGLState&, void> updateColorFunc);
  
  virtual void initGTK();

  virtual void showControls(Gtk::ScrolledWindow* win);
    
  typedef enum {
    SINGLE_COLOR = 1,
    COLOR_BY_ID = 2,
    COLOR_BY_SPEED = 3
  } DrawMode;

  inline volatile const DrawMode& getDrawMode() { return _mode; }
  inline volatile const double& getScaleV() { return _scaleV; }
  inline volatile const cl_uchar4& getColorFixed() { return _colorFixed; }
  inline volatile const cl_uchar4& getColorStatic() { return _colorStatic; }
  inline volatile const bool& getRecolorOnUpdate() { return _recolorOnUpdate; }
  inline volatile bool getColorIfStatic() { return _colorStaticParticles; }
protected:

  void guiUpdate();

  std::auto_ptr<Gtk::VBox> _optList;
  std::auto_ptr<Gtk::CheckButton> _colorIfStatic;
  std::auto_ptr<Gtk::SpinButton> _RStatic;
  std::auto_ptr<Gtk::SpinButton> _GStatic;
  std::auto_ptr<Gtk::SpinButton> _BStatic;
  std::auto_ptr<Gtk::SpinButton> _AStatic;
  std::auto_ptr<Gtk::RadioButton> _singleColorMode;
  std::auto_ptr<Gtk::RadioButton> _colorByIDMode;
  std::auto_ptr<Gtk::RadioButton> _colorBySpeedMode;
  std::auto_ptr<Gtk::SpinButton> _RFixed;
  std::auto_ptr<Gtk::SpinButton> _GFixed;
  std::auto_ptr<Gtk::SpinButton> _BFixed;
  std::auto_ptr<Gtk::SpinButton> _AFixed;
  std::auto_ptr<Gtk::Entry> _characteristicV;
    
  volatile cl_uchar4 _colorFixed;
  volatile cl_uchar4 _colorStatic;
  volatile DrawMode _mode;
  volatile double _scaleV;
  volatile bool _recolorOnUpdate;
  volatile bool _colorStaticParticles;

  magnet::function::Delegate1<magnet::CL::CLGLState&, void> _updateColorFunc;
};
#endif

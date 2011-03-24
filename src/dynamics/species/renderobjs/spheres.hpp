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
# include <magnet/gtk/colorMapSelector.hpp>

class SphereParticleRenderer: public RTSpheres
{
public:
  //! \param spheresPerObject Used when one simulation object is
  //! represented by many spheres. If r_{i,a} is the a'th sphere of
  //! object i, then _particleData should contain coordinates like
  //! (r_(0,0),r_(1,0),r_(2,0)....r_(0,1),r_(1,1)...). But only
  //! N/_spheresPerObject colors should be placed in
  //! _particleColorData and the data is duplicated out automatically
  //! to all spheres in a single object.
  SphereParticleRenderer(size_t N, std::string name, 
			 magnet::function::Delegate1<magnet::CL::CLGLState&, void> updateColorFunc,
			 size_t spheresPerObject = 1);
  
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
  inline volatile bool getRecolorOnUpdate() { return _colorStaticParticles || _recolorOnUpdate; }
  inline volatile bool getColorIfStatic() { return _colorStaticParticles; }

  inline void map(cl_uchar4 &color, float val) { _colorMap->map(color, val); }

  std::vector<cl_float4> _particleData;
  std::vector<cl_uchar4> _particleColorData;

  void sendRenderData(magnet::CL::CLGLState& CLState);
  void sendColorData(magnet::CL::CLGLState& CLState);

protected:
  size_t _spheresPerObject;

  void guiUpdate();

  std::auto_ptr<Gtk::VBox> _optList;
  std::auto_ptr<magnet::Gtk::ColorMapSelector> _colorMap;
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

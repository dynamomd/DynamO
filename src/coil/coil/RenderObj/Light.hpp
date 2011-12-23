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

#include "RenderObj.hpp"
#include <magnet/gtk/transferFunction.hpp>
#include <magnet/GL/texture.hpp>
#include <magnet/GL/shader/volume.hpp>
#include <magnet/GL/objects/cube.hpp>
#include <memory>
#include <tr1/array>

namespace coil {
  class RLight : public RenderObj
  {
  public:
    RLight(std::string name): 
      RenderObj(name), _intensity(1), 
      _attenuation(0.5), _specularExponent(96),
      _specularFactor(1) {}
  
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();
    virtual void forwardRender(magnet::GL::FBO& fbo,
			       const magnet::GL::Camera& cam,
			       const magnet::GL::Light& light,
			       RenderMode mode) {}

    virtual void clTick(const magnet::GL::Camera&) {}

    virtual void showControls(Gtk::ScrolledWindow* win);

    float getIntensity() const { return _intensity; }
    float getAttenuation() const { return _attenuation; }
    float getSpecularExponent() const { return _specularExponent; }
    float getSpecularFactor() const { return _specularFactor; }

  protected:
    void initGTK();
    void guiUpdate();

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    //GTK gui stuff
    std::auto_ptr<Gtk::VBox> _optList;
    std::auto_ptr<Gtk::Entry> _intensityEntry;
    std::auto_ptr<Gtk::Entry> _attenuationEntry;
    std::auto_ptr<Gtk::Entry> _specularExponentEntry;
    std::auto_ptr<Gtk::Entry> _specularFactorEntry;
    
    float _intensity, _attenuation, _specularExponent, _specularFactor;
  };
}

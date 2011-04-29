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
#include <magnet/GL/SSAO.hpp>
#include "filter.hpp"

namespace coil 
{
  class SSAOWrapper: public filter
  {
  public:
    SSAOWrapper();
    ~SSAOWrapper();

    inline virtual size_t type_id() { return detail::filterEnum<SSAOWrapper>::val; }    

    inline virtual void invoke(GLuint colorTextureUnit, size_t width, size_t height,
			       const magnet::GL::viewPort& vp);

    inline virtual bool needsNormalDepth()  { return true; }

    virtual void showControls(Gtk::ScrolledWindow*);

  protected:
    magnet::GL::SSAO _filter;
    Glib::RefPtr<Gtk::Builder> _refXml;
    GLfloat _radius;
    GLfloat _totStrength;
    GLfloat _dropoff;

    GLuint _randomTexture;

    void settingsCallback();
    
    Gtk::HScale _radiusSlider;
    Gtk::HScale _totStrengthSlider;
    Gtk::HScale _dropoffSlider;
    Gtk::HBox _optlist;
    static const size_t _randomTextureSize = 64;
  };
}

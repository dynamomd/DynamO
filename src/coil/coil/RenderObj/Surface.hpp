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

#include "Triangles.hpp"
#include <time.h>
#include <magnet/math/vector.hpp>

namespace coil {
  class RSurface : public RTriangles
  {
  public:
    RSurface(magnet::GL::Context::ContextPtr context, std::string name, size_t N = 10, Vector origin = Vector{-25,-1.5,-25}, Vector axis1 = Vector{50,0,0},
	     Vector axis2 = Vector{0,0,50}, Vector axis3 = Vector{0,1,0});

    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    virtual bool deletable() { return true; }

    virtual magnet::math::Vector getMaxCoord() const;
    virtual magnet::math::Vector getMinCoord() const;

  protected:
    void clTick() {}

    size_t _N;

    Vector _origin;
    Vector _axis1;
    Vector _axis2;     
    Vector _axis3;
  };
}

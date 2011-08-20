/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include <dynamo/dynamics/globals/global.hpp>
#include <vector>
#include <map>

namespace dynamo {
  class GWaker: public Global
  {
  public:
    GWaker(const magnet::xml::Node&, dynamo::SimData*);

    GWaker(dynamo::SimData*, const std::string&, CRange*, const double, const double,
	   std::string nblist);
  
    virtual ~GWaker() {}

    virtual Global* Clone() const { return new GWaker(*this); };

    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(const Particle&, const double) const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

  protected:
    void particlesUpdated(const NEventData&);

    void nblistCallback(const Particle& part, const size_t& oid) const;

    mutable size_t _neighbors;

    virtual void outputXML(magnet::xml::XmlStream&) const;
    double _wakeTime;
    double _wakeVelocity;

    std::string _nblistName;
    size_t _NBListID;  
  };
}

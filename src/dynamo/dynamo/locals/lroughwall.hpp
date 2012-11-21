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
#include <dynamo/locals/local.hpp>

namespace dynamo {
  class LRoughWall: public Local
  {
  public:
    LRoughWall(const magnet::xml::Node&, dynamo::Simulation*);
    LRoughWall(dynamo::Simulation*, double, double, double, Vector , Vector , 
	       std::string, IDRange*, bool nrender = true);

    virtual ~LRoughWall() {}

    virtual LocalEvent getEvent(const Particle&) const;

    virtual void runEvent(Particle&, const LocalEvent&) const;
  
    virtual void operator<<(const magnet::xml::Node&);

    virtual void checkOverlaps(const Particle&) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    Vector  vNorm;
    Vector  vPosition;
    double e;
    double et;
    double r;
    bool render;
  };
}

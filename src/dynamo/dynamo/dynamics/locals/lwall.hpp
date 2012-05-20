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
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/base/is_simdata.hpp>

namespace dynamo {
  class LWall: public Local
  {
  public:
    LWall(const magnet::xml::Node&, dynamo::SimData*);

    template<class T1, class T2>
    LWall(dynamo::SimData* nSim, T1 e, T2 d, Vector nnorm,
	  Vector  norigin, std::string nname, Range* nRange):
      Local(nRange, nSim, "LocalWall"),
      vNorm(nnorm),
      vPosition(norigin),
      _diameter(Sim->_properties.getProperty
		(d, Property::Units::Length())),
      _e(Sim->_properties.getProperty
	 (e, Property::Units::Dimensionless()))
    { localName = nname; }

    virtual ~LWall() {}

    virtual LocalEvent getEvent(const Particle&) const;

    virtual void runEvent(const Particle&, const LocalEvent&) const;
  
    virtual bool isInCell(const Vector &, const Vector &) const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

    virtual void checkOverlaps(const Particle&) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    Vector  vNorm;
    Vector  vPosition;
    shared_ptr<Property> _diameter;
    shared_ptr<Property> _e;
    bool render;
  };
}

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
#include "local.hpp"
#include <boost/tuple/tuple.hpp>
#include <vector>

class LTriangleMesh: public Local
{
public:
  LTriangleMesh(const magnet::xml::Node&, dynamo::SimData*);
  LTriangleMesh(dynamo::SimData*, double e, std::string, CRange*);

  virtual ~LTriangleMesh() {}

  virtual Local* Clone() const { return new LTriangleMesh(*this); };

  virtual LocalEvent getEvent(const Particle&) const;

  virtual void runEvent(const Particle&, const LocalEvent&) const;
  
  virtual bool isInCell(const Vector &, const Vector &) const;

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node&);

  virtual void checkOverlaps(const Particle&) const;

protected:
  virtual void outputXML(xml::XmlStream&) const;

  std::vector<Vector> _vertices;

  typedef boost::tuples::tuple<size_t, size_t, size_t> TriangleElements;
  std::vector<TriangleElements> _elements;

  double _e;
};

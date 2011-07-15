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
#include "../coilRenderObj.hpp"
#include "../../base/is_simdata.hpp"
#include <boost/tuple/tuple.hpp>
#include <vector>


class LTriangleMesh: public Local, public CoilRenderObj
{
public:
  LTriangleMesh(const magnet::xml::Node&, dynamo::SimData*);

  template<class T1, class T2>
  LTriangleMesh(dynamo::SimData* nSim, T1 e, T2 d, std::string name, CRange* nRange):
    Local(nRange, nSim, "LocalWall"),
    _e(Sim->_properties.getProperty
       (e, Property::Units::Dimensionless())),
    _diameter(Sim->_properties.getProperty
	      (d, Property::Units::Length()))
  { localName = name; }

  virtual ~LTriangleMesh() {}

  virtual Local* Clone() const { return new LTriangleMesh(*this); };

  virtual LocalEvent getEvent(const Particle&) const;

  virtual void runEvent(const Particle&, const LocalEvent&) const;
  
  virtual bool isInCell(const Vector &, const Vector &) const;

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node&);

  virtual void checkOverlaps(const Particle&) const;

#ifdef DYNAMO_visualizer
  virtual magnet::thread::RefPtr<RenderObj>& getCoilRenderObj() const;
  virtual void updateRenderData(magnet::GL::Context&) const;
#endif

protected:
#ifdef DYNAMO_visualizer
  mutable magnet::thread::RefPtr<RenderObj> _renderObj;
#endif

  virtual void outputXML(magnet::xml::XmlStream&) const;

  std::vector<Vector> _vertices;

  typedef boost::tuples::tuple<size_t, size_t, size_t> TriangleElements;
  std::vector<TriangleElements> _elements;

  magnet::thread::RefPtr<Property> _e;
  magnet::thread::RefPtr<Property> _diameter;
};

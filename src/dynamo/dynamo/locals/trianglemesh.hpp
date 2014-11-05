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
#include <dynamo/coilRenderObj.hpp>
#include <dynamo/simulation.hpp>
#include <tuple>
#include <vector>

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/TriangleMesh.hpp>
#endif

namespace dynamo {
  class LTriangleMesh: public Local, public CoilRenderObj
  {
  public:
    LTriangleMesh(const magnet::xml::Node&, dynamo::Simulation*);

    template<class T1, class T2>
    LTriangleMesh(dynamo::Simulation* nSim, T1 e, T2 d, std::string name, IDRange* nRange):
      Local(nRange, nSim, "LocalWall"),
      _e(Sim->_properties.getProperty(e, Property::Units::Dimensionless())),
      _diameter(Sim->_properties.getProperty(d, Property::Units::Length()))
    { localName = name; }

    virtual ~LTriangleMesh() {}

    virtual Event getEvent(const Particle&) const;

    virtual ParticleEventData runEvent(Particle&, const Event&) const;
  
    virtual void operator<<(const magnet::xml::Node&);

    virtual bool validateState(const Particle& part, bool textoutput = true) const { return false; }

#ifdef DYNAMO_visualizer
    virtual shared_ptr<coil::RenderObj> getCoilRenderObj() const;
    virtual void updateRenderData() const {}
#endif

  protected:
#ifdef DYNAMO_visualizer
    mutable shared_ptr<coil::RTriangleMesh> _renderObj;
#endif

    virtual void outputXML(magnet::xml::XmlStream&) const;

    std::vector<Vector> _vertices;

    typedef std::tuple<size_t, size_t, size_t> TriangleElements;
    std::vector<TriangleElements> _elements;

    shared_ptr<Property> _e;
    shared_ptr<Property> _diameter;
  };
}

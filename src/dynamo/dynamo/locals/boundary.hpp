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
#include <dynamo/locals/boundary/object.hpp>
#include <dynamo/coilRenderObj.hpp>
#ifdef DYNAMO_visualizer
# include <coil/RenderObj/TriangleMesh.hpp>
#endif 

namespace dynamo {
  class LBoundary: public Local, public CoilRenderObj
  {
  public:
    LBoundary(const magnet::xml::Node&, dynamo::Simulation*);

    virtual ~LBoundary() {}

    virtual Event getEvent(const Particle&) const;

    virtual ParticleEventData runEvent(Particle&, const Event&) const;
  
    virtual void operator<<(const magnet::xml::Node&);

    virtual bool validateState(const Particle& part, bool textoutput = true) const;

#ifdef DYNAMO_visualizer
    std::pair<std::vector<float>, std::vector<GLuint> > getTessalatedSurfaces() const;
    virtual shared_ptr<coil::RenderObj> getCoilRenderObj() const;
    virtual void updateRenderData() const;
#endif

  protected:
#ifdef DYNAMO_visualizer
    mutable shared_ptr<coil::RTriangleMesh> _renderObj;
#endif

    virtual void outputXML(magnet::xml::XmlStream&) const;

    shared_ptr<Property> _diameter;
    double _kT;
    std::vector<std::shared_ptr<boundary::Object> > _objects;

    boundary::BoundaryOscillationData _oscillationData;
  };
}

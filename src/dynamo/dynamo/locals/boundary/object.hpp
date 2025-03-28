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
#include <dynamo/1particleEventData.hpp>
#include <dynamo/base.hpp>
#include <magnet/math/vector.hpp>
#ifdef DYNAMO_visualizer
#include <coil/RenderObj/TriangleMesh.hpp>
#endif

namespace magnet {
namespace xml {
class XmlStream;
}
} // namespace magnet

namespace dynamo {
class Particle;
class Event;

namespace boundary {

struct BoundaryOscillationData {
  Vector _origin;
  Vector _amplitude;
  double _freq;
  double _t_shift;
};

class Object : public dynamo::SimBase {
public:
  virtual ~Object() {}

  Object(Simulation *const SD, const std::string aName,
         const BoundaryOscillationData &data);

  virtual bool validateState(const Particle &part, bool textoutput) const = 0;

  virtual Event getEvent(const Particle &part, const double diameter) const = 0;

  virtual Vector getContactNormal(const Particle &, const Event &) const = 0;

  virtual void outputXML(magnet::xml::XmlStream &) const = 0;

  static shared_ptr<Object> getClass(const magnet::xml::Node &,
                                     dynamo::Simulation *,
                                     const BoundaryOscillationData &);

#ifdef DYNAMO_visualizer
  virtual std::pair<std::vector<float>, std::vector<GLuint>>
  getTessalatedSurfaces() const = 0;
#endif
protected:
  const BoundaryOscillationData &_oscillationData;
};

class PlanarWall : public Object {
public:
  PlanarWall(const magnet::xml::Node &, dynamo::Simulation *,
             const BoundaryOscillationData &data);

  virtual bool validateState(const Particle &part,
                             bool textoutput = true) const;

  virtual Event getEvent(const Particle &part, const double diameter) const;

  virtual Vector getContactNormal(const Particle &, const Event &) const;

  virtual void outputXML(magnet::xml::XmlStream &) const;

#ifdef DYNAMO_visualizer
  std::pair<std::vector<float>, std::vector<GLuint>>
  getTessalatedSurfaces() const;
#endif

protected:
  Vector _position;
  Vector _normal;
};
} // namespace boundary
} // namespace dynamo

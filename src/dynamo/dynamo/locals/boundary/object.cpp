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
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/locals/boundary/object.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
namespace boundary {
Object::Object(Simulation *const SD, const std::string aName,
               const BoundaryOscillationData &data)
    : dynamo::SimBase(SD, aName), _oscillationData(data) {}

shared_ptr<Object> Object::getClass(const magnet::xml::Node &XML,
                                    dynamo::Simulation *Sim,
                                    const BoundaryOscillationData &data) {
  if (!XML.getAttribute("Type").getValue().compare("PlanarWall"))
    return shared_ptr<Object>(new PlanarWall(XML, Sim, data));
  else
    M_throw() << XML.getAttribute("Type").getValue()
              << ", Unknown type of Object encountered" << XML.getPath();
}

PlanarWall::PlanarWall(const magnet::xml::Node &XML, dynamo::Simulation *Sim,
                       const BoundaryOscillationData &data)
    : Object(Sim, "PlanarWall", data) {
  if (XML.hasNode("Position"))
    _position << XML.getNode("Position");
  _position *= Sim->units.unitLength();
  _normal << XML.getNode("Normal");
  _normal.normalise();
}

void PlanarWall::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "PlanarWall";
  if (_position.nrm() != 0)
    XML << magnet::xml::tag("Position") << _position / Sim->units.unitLength()
        << magnet::xml::endtag("Position");

  XML << magnet::xml::tag("Normal") << _normal << magnet::xml::endtag("Normal");
}

bool PlanarWall::validateState(const Particle &part, bool textoutput) const {
  return true;
}

Event PlanarWall::getEvent(const Particle &part, const double diameter) const {
  const double dt = Sim->dynamics->getPlaneEvent(
      part, _position + _oscillationData._origin, _normal, diameter / 2);
  return Event(part.getID(), dt, LOCAL,
               (dt == std::numeric_limits<float>::infinity()) ? NONE : WALL, 0);
}

Vector PlanarWall::getContactNormal(const Particle &, const Event &) const {
  return _normal;
}

#ifdef DYNAMO_visualizer
std::pair<std::vector<float>, std::vector<GLuint>>
PlanarWall::getTessalatedSurfaces() const {
  // Intersect the plane of the surface, with the unit cell box to
  // generate a polygon.

  // This approach is detailed in "A vertex program for efficient
  // box-plane intersection" by Christof Rezk Salama and Andreas
  // Kolb.

  // Generate the vertices of the cube, in a certain ordering.

  // Each vertex can be represented by three bools. Each bool
  // indicates if the position is at the minimum (false) or the
  // maximum (true) in that direction. This is assuming an axis aligned cube.
  std::array<bool, 3> frontvertex = {{false, false, false}};

  auto vertexPos = [&](const std::array<bool, 3> &vertex) {
    return magnet::math::elementwiseMultiply(
        Vector{float(vertex[0]), float(vertex[1]), float(vertex[2])} -
            Vector{0.5, 0.5, 0.5},
        Sim->primaryCellSize);
  };

  // We use the max distance to determine the vertex which is the
  // furthest in front of the plane. We use the min_distance to
  // check that the cube and plane intersect.
  double max_distance = -std::numeric_limits<float>::infinity();
  double min_distance = +std::numeric_limits<float>::infinity();
  for (size_t i(0); i < 8; ++i) {
    std::array<bool, 3> vertex = {
        {bool((i >> 0) & 1), bool((i >> 1) & 1), bool((i >> 2) & 1)}};

    const Vector vertexpos = vertexPos(vertex);

    const double vertexdistance = (vertexpos - _position) * _normal;
    if (vertexdistance > max_distance) {
      max_distance = vertexdistance;
      frontvertex = vertex;
    }
    min_distance = std::min(vertexdistance, min_distance);
  }

  if ((min_distance > 0) || (max_distance < 0))
    M_throw() << "Cannot correctly render a wall which lies outside of the "
                 "primary image!";

  // We genertate an ordering of vertices corresponding to the
  // example in Fig. 3.
  std::array<std::array<bool, 3>, 8> vertex_bits = {
      {{{
           frontvertex[0],
           frontvertex[0],
           frontvertex[0],
       }},
       {{!frontvertex[0], frontvertex[0], frontvertex[0]}},
       {{frontvertex[0], !frontvertex[0], frontvertex[0]}},
       {{frontvertex[0], frontvertex[0], !frontvertex[0]}},
       {{!frontvertex[0], frontvertex[0], !frontvertex[0]}},
       {{!frontvertex[0], !frontvertex[0], frontvertex[0]}},
       {{frontvertex[0], !frontvertex[0], !frontvertex[0]}},
       {{!frontvertex[0], !frontvertex[0], !frontvertex[0]}}}};

  std::array<Vector, 8> vertex_pos;
  for (size_t i(0); i < 8; ++i)
    vertex_pos[i] = vertexPos(vertex_bits[i]);

  const double d = _normal * (_position + _oscillationData._origin);

  auto lambda = [&](const size_t id1, const size_t id2) {
    return (d - _normal * vertex_pos[id1]) /
           (_normal * (vertex_pos[id2] - vertex_pos[id1]));
  };

  Vector P0;
  for (auto id_pair :
       std::array<std::array<size_t, 2>, 3>{{{{0, 1}}, {{1, 4}}, {{4, 7}}}}) {
    const double l = lambda(id_pair[0], id_pair[1]);
    if ((l >= 0) && (l <= 1)) {
      P0 = (1 - l) * vertex_pos[id_pair[0]] + l * vertex_pos[id_pair[1]];
      break;
    }
  }

  Vector P2;
  for (auto id_pair :
       std::array<std::array<size_t, 2>, 3>{{{{0, 2}}, {{2, 5}}, {{5, 7}}}}) {
    double l = lambda(id_pair[0], id_pair[1]);
    if ((l >= 0) && (l <= 1)) {
      P2 = (1 - l) * vertex_pos[id_pair[0]] + l * vertex_pos[id_pair[1]];
      break;
    }
  }

  Vector P4;
  for (auto id_pair :
       std::array<std::array<size_t, 2>, 3>{{{{0, 3}}, {{3, 6}}, {{6, 7}}}}) {
    double l = lambda(id_pair[0], id_pair[1]);
    if ((l >= 0) && (l <= 1)) {
      P4 = (1 - l) * vertex_pos[id_pair[0]] + l * vertex_pos[id_pair[1]];
      break;
    }
  }

  Vector P1 = P0;
  {
    const double l = lambda(1, 5);
    if ((l >= 0) && (l <= 1)) {
      P1 = (1 - l) * vertex_pos[1] + l * vertex_pos[5];
    }
  }

  Vector P3 = P2;
  {
    const double l = lambda(2, 6);
    if ((l >= 0) && (l <= 1)) {
      P3 = (1 - l) * vertex_pos[2] + l * vertex_pos[6];
    }
  }

  Vector P5 = P4;
  {
    const double l = lambda(3, 4);
    if ((l >= 0) && (l <= 1)) {
      P5 = (1 - l) * vertex_pos[3] + l * vertex_pos[4];
    }
  }

  const std::vector<float> vertices = {
      float(P0[0]), float(P0[1]), float(P0[2]), float(P1[0]), float(P1[1]),
      float(P1[2]), float(P2[0]), float(P2[1]), float(P2[2]), float(P3[0]),
      float(P3[1]), float(P3[2]), float(P4[0]), float(P4[1]), float(P4[2]),
      float(P5[0]), float(P5[1]), float(P5[2])};

  const std::vector<GLuint> indices = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5};
  return std::make_pair(vertices, indices);
}
#endif
} // namespace boundary
} // namespace dynamo

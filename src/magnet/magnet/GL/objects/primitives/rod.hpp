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
#include <cmath>
#include <vector>

#include <GL/gl.h>
#include <magnet/exception.hpp>
#include <magnet/math/vector.hpp>
#include <stdlib.h>
#include <string.h>

namespace magnet {
namespace GL {
namespace objects {
namespace primitives {
/*! \brief This class contains functions which generate the
    vertex data for an OpenGL rod (cylinder with closed ends).
 */
class Rod {
public:
  inline static std::vector<GLfloat> getVertices(size_t LOD) {
    std::vector<GLfloat> vertices;
    vertices.reserve((6 * LOD + 2) * 3);

    for (size_t vert = 0; vert < 2 * LOD; ++vert) {
      vertices.push_back(0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD));
      vertices.push_back(0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD));
      vertices.push_back(vert % 2 - 0.5);
    }

    // Put vertices around the edge of each side of the
    // cylinder (we need to repeat these from the ones above to
    // get the normals correct for the end faces)
    for (size_t vert = 0; vert < LOD; ++vert) {
      vertices.push_back(0.5f * std::sin(vert * 2.0f * M_PI / LOD));
      vertices.push_back(0.5f * std::cos(vert * 2.0f * M_PI / LOD));
      vertices.push_back(-0.5f);
    }

    vertices.push_back(0);
    vertices.push_back(0);
    vertices.push_back(-0.5f);

    for (size_t vert = 0; vert < LOD; ++vert) {
      vertices.push_back(0.5f * std::sin(vert * 2.0f * M_PI / LOD));
      vertices.push_back(0.5f * std::cos(vert * 2.0f * M_PI / LOD));
      vertices.push_back(0.5f);
    }

    vertices.push_back(0);
    vertices.push_back(0);
    vertices.push_back(0.5f);

    return vertices;
  }

  inline static std::vector<GLfloat> getNormals(size_t LOD) {
    std::vector<GLfloat> normals;
    normals.reserve((6 * LOD + 2) * 3);

    for (size_t vert = 0; vert < 2 * LOD; ++vert) {
      GLfloat x = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD);
      GLfloat y = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD);
      GLfloat scale = 1.0f / std::sqrt(x * x + y * y);
      normals.push_back(x * scale);
      normals.push_back(y * scale);
      normals.push_back(0);
    }

    for (size_t vert = 0; vert < LOD + 1; ++vert) {
      normals.push_back(0);
      normals.push_back(0);
      normals.push_back(-1.0);
    }

    for (size_t vert = 0; vert < LOD + 1; ++vert) {
      normals.push_back(0);
      normals.push_back(0);
      normals.push_back(1.0);
    }

    return normals;
  }

  inline static std::vector<GLuint> getIndices(size_t LOD) {
    // 2 triangles per face (3 indices per triangle)
    std::vector<GLuint> indices;
    for (size_t vert = 0; vert < LOD; ++vert) {
      indices.push_back((2 * vert + 0) % (2 * LOD));
      indices.push_back((2 * vert + 1) % (2 * LOD));
      indices.push_back((2 * vert + 2) % (2 * LOD));
      indices.push_back((2 * vert + 1) % (2 * LOD));
      indices.push_back((2 * vert + 3) % (2 * LOD));
      indices.push_back((2 * vert + 2) % (2 * LOD));
    }

    for (size_t vert = 0; vert < LOD; ++vert) {
      indices.push_back(2 * LOD + vert);
      indices.push_back(2 * LOD + ((vert + 1) % LOD));
      indices.push_back(3 * LOD + 0);

      indices.push_back(3 * LOD + 1 + ((vert + 1) % LOD));
      indices.push_back(3 * LOD + 1 + vert);
      indices.push_back(4 * LOD + 1);
    }

    // End faces
    return indices;
  }

protected:
};
} // namespace primitives
} // namespace objects
} // namespace GL
} // namespace magnet

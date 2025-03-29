/*    dynamo:- Event driven molecular dynamics simulator
 *    http://www.dynamomd.org
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#include <magnet/GL/shader/detail/filter.hpp>

namespace magnet {
namespace GL {
namespace shader {
/*! \brief Implements a 5x5 Laplacian filter Shader. */
class Laplacian5x5 : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(5); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[5][5] = {{0, 0, -1, 0, 0},
                                          {0, -1, -2, -1, 0},
                                          {-1, -2, 16, -2, -1},
                                          {0, -1, -2, -1, 0},
                                          {0, 0, -1, 0, 0}};

    return (const GLfloat *)weights;
  }
};

/*! \brief Implements variant A of the 3x3 Laplacian filter Shader. */
class Laplacian3x3A : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(3); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[3][3] = {{0, -1, 0}, {-1, 4, -1}, {0, -1, 0}};

    return (const GLfloat *)weights;
  }
};

/*! \brief Implements variant B of the 3x3 Laplacian filter Shader. */
class Laplacian3x3B : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(3); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[3][3] = {
        {-1, -1, -1}, {-1, 8, -1}, {-1, -1, -1}};

    return (const GLfloat *)weights;
  }
};

/*! \brief Implements a 9x9 Laplacian of Gaussian Blurred Shader. */
class LoG9x9 : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(9); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[9][9] = {
        {0, 1, 1, 2, 2, 2, 1, 1, 0},       {1, 2, 4, 5, 5, 5, 4, 2, 1},
        {1, 4, 5, 3, 0, 3, 5, 4, 1},       {2, 5, 3, -12, -24, -12, 3, 5, 2},
        {2, 5, 0, -24, -40, -24, 0, 5, 2}, {2, 5, 3, -12, -24, -12, 3, 5, 2},
        {1, 4, 5, 3, 0, 3, 5, 4, 1},       {1, 2, 4, 5, 5, 5, 4, 2, 1},
        {0, 1, 1, 2, 2, 2, 1, 1, 0}};

    return (const GLfloat *)weights;
  }
};
} // namespace shader
} // namespace GL
} // namespace magnet

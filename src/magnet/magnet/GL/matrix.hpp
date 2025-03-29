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
#include <array>
#include <magnet/math/matrix.hpp>

namespace magnet {
namespace GL {
typedef magnet::math::NMatrix<GLfloat, 4> GLMatrix;

inline math::Matrix demoteToMatrix(const GLMatrix &M) {
  math::Matrix retval;
  for (size_t i(0); i < 3; ++i)
    for (size_t j(0); j < 3; ++j)
      retval(i, j) = M(i, j);
  return retval;
}

inline GLMatrix promoteToGLMatrix(const math::Matrix &M) {
  GLMatrix retval;
  for (size_t i(0); i < 3; ++i)
    for (size_t j(0); j < 3; ++j)
      retval(i, j) = M(i, j);
  retval(3, 3) = 1;
  return retval;
}

/*! \brief Return a matrix corresponding to a translation.
 *
 * This command emulates the glTranslate command.
 */
inline GLMatrix translate(GLfloat x, GLfloat y, GLfloat z) {
  return GLMatrix{1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1};
}
inline GLMatrix translate(const math::Vector &vec) {
  return translate(vec[0], vec[1], vec[2]);
}

/*! \brief Return a matrix corresponding to a scaling.
 *
 * This command emulates the glScale command.
 */
inline GLMatrix scale(GLfloat x, GLfloat y, GLfloat z) {
  return GLMatrix{x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1};
}
/*! \brief Return a matrix corresponding to a scaling.
 *
 * This command emulates the glScale command.
 */
inline GLMatrix scale(const math::Vector &vec) {
  return scale(vec[0], vec[1], vec[2]);
}

/*! \brief Return a matrix corresponding to a rotation.
 *
 * This command emulates the glRotate command.
 *
 * \param angle The angle of rotation (in degrees).
 * \param axis The axis of rotation.
 */
inline GLMatrix rotate(const GLfloat &angle, const math::Vector &axis) {
  return promoteToGLMatrix(Rodrigues((angle * M_PI / 180.0f) * axis));
}

/*! \brief Return a matrix corresponding to a frustrum projection.

  This command emulates the glFrustrum command with one
  important exception. There is an additional factor called
  zoffset, which biases all surfaces towards (positive) or away
  (negative) from the camera. This is used to solve Z-fighting
  errors. The resource which explains this value is given here
  http://www.terathon.com/gdc07_lengyel.pdf

  If you wish to bias a light source's projection matrix (for
  shadow map calculations) you should set zoffset to 4.8e-7.
*/
inline GLMatrix frustrum(const GLfloat left, const GLfloat right,
                         const GLfloat bottom, const GLfloat top,
                         const GLfloat nearVal, const GLfloat farVal,
                         const GLfloat zoffset = 0) {
  GLfloat A = (right + left) / (right - left);
  GLfloat B = (top + bottom) / (top - bottom);
  GLfloat C = -(farVal + nearVal) / (farVal - nearVal);
  GLfloat D = -2 * farVal * nearVal / (farVal - nearVal);

  return GLMatrix{2 * nearVal / (right - left),
                  0,
                  A,
                  0,
                  0,
                  2 * nearVal / (top - bottom),
                  B,
                  0,
                  0,
                  0,
                  C,
                  D,
                  0,
                  0,
                  -1,
                  0};
}
/*! \brief Return a matrix corresponding to a perspective projection.
 *
 * This command emulates the gluPerspective command.
 */
inline GLMatrix perspective(const GLfloat fovy, const GLfloat aspect,
                            const GLfloat zNear, const GLfloat zFar) {
  GLfloat f = 1 / std::tan(fovy * 0.5f);
  return GLMatrix{f / aspect,
                  0,
                  0,
                  0,
                  0,
                  f,
                  0,
                  0,
                  0,
                  0,
                  (zFar + zNear) / (zNear - zFar),
                  2 * zFar * zNear / (zNear - zFar),
                  0,
                  0,
                  -1,
                  0};
}
} // namespace GL
} // namespace magnet

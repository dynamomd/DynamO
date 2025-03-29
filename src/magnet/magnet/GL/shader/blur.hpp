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
#define STRINGIFY(A) #A

namespace magnet {
namespace GL {
namespace shader {

/*! \brief A seperable Gaussian blur equilivent to a 13x13 kernel.

  This Gaussian kernel is discussed here
  http://rastergrid.com/blog/2010/09/efficient-gaussian-blur-with-linear-sampling/

  The input texture must have linear sampling enabled.
 */
class SeperableGaussian : public detail::SSShader {
public:
  virtual std::string initFragmentShaderSource() {
    return STRINGIFY(
        layout(location = 0) out vec4 color_out;

        // The HDR color buffer
        uniform sampler2D colorTex; uniform vec2 invDim; uniform int direction;

        const float offset[3] = float[](0.0, 1.3846153846, 3.2307692308);
        const float weight[3] =
            float[](0.2270270270, 0.3162162162, 0.0702702703);

        void main() {
          vec4 sum = weight[0] * texture(colorTex, gl_FragCoord.xy * invDim);

          if (direction == 1)
            for (int i = 1; i < 3; ++i)
              sum +=
                  weight[i] *
                  (texture(colorTex,
                           (gl_FragCoord.xy + vec2(0.0, offset[i])) * invDim) +
                   texture(colorTex,
                           (gl_FragCoord.xy - vec2(0.0, offset[i])) * invDim));
          else
            for (int i = 1; i < 3; ++i)
              sum +=
                  weight[i] *
                  (texture(colorTex,
                           (gl_FragCoord.xy + vec2(offset[i], 0.0)) * invDim) +
                   texture(colorTex,
                           (gl_FragCoord.xy - vec2(offset[i], 0.0)) * invDim));

          color_out = sum;
        });
  }
};

/*! \brief Implements a 5x5 Gaussian Blur Shader.
 */
class Gaussian5x5Blur : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(5); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[5][5] = {
        {1.0 / 331.0, 4 / 331.0, 7 / 331.0, 4 / 331.0, 1 / 331.0},
        {4 / 331.0, 20 / 331.0, 33 / 331.0, 20 / 331.0, 4 / 331.0},
        {7 / 331.0, 33 / 331.0, 55 / 331.0, 33 / 331.0, 7 / 331.0},
        {4 / 331.0, 20 / 331.0, 33 / 331.0, 20 / 331.0, 4 / 331.0},
        {1 / 331.0, 4 / 331.0, 7 / 331.0, 4 / 331.0, 1 / 331.0}};

    return (const GLfloat *)weights;
  }
};

/*! \brief Implements a 5x5 Box Blur Shader. */
class Box5x5Blur : public detail::SSKernelShader {
public:
  void build() { SSKernelShader::build(5); }

  virtual const GLfloat *weights() {
    static const GLfloat weights[5][5] = {
        {1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0},
        {1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0},
        {1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0},
        {1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0},
        {1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0, 1 / 25.0}};

    return (const GLfloat *)weights;
  }
};
} // namespace shader
} // namespace GL
} // namespace magnet

#undef STRINGIFY

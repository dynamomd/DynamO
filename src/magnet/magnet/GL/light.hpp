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

//Here we have the correct order of GL includes
#include <magnet/GL/context.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/actor.hpp>

namespace magnet
{
    namespace GL
    {
        class Light : public Actor, public magnet::GL::CameraHeadTracking
        {
        protected:
            bool _shadow_casting = false;
            float _intensity, _specularExponent, _specularFactor;
            float _maxvariance, _bleedreduction;
            GLfloat _size;
            std::array<GLfloat, 3> _color;

        public:
            Light(magnet::GL::Context::ContextPtr context, magnet::math::Vector position, magnet::math::Vector lookAtPoint,
                  GLfloat zNearDist = 8.0f, GLfloat zFarDist = 10000.0f, magnet::math::Vector up = magnet::math::Vector{0, 1, 0},
                  GLfloat simLength = 25.0f, GLfloat size = 0.2) : CameraHeadTracking(context, position, lookAtPoint, zNearDist, zFarDist, up, simLength, magnet::math::Vector{0, 0, 20}),
                                                                   _shadow_casting(false),
                                                                   _intensity(1.0),
                                                                   _specularExponent(32),
                                                                   _specularFactor(1),
                                                                   _maxvariance(0.1),
                                                                   _bleedreduction(0.2),
                                                                   _size(size),
                                                                   _color{1.0, 1.0, 1.0}
            {
            }

            virtual ~Light() {}

            float getIntensity() const { return _intensity; }
            float getSpecularExponent() const { return _specularExponent; }
            float getSpecularFactor() const { return _specularFactor; }
            float getMaxVariance() const { return _maxvariance; }
            float getBleedReduction() const { return _bleedreduction; }
            bool shadowCasting() const { return _shadow_casting; }

            /*! \brief Load the specified OpenGL texture matrix with the
      projection required for shadow mapping.
      
      \note The current OpenGL model view matrix must be the matrix
      used for rendering.
    */
            inline magnet::GL::GLMatrix getShadowTextureMatrix()
            {
                return magnet::GL::translate(magnet::math::Vector{0.5, 0.5, 0.5}) * magnet::GL::scale(magnet::math::Vector{0.5, 0.5, 0.5}) * getProjectionMatrix() * getViewMatrix();
            }

            const std::array<GLfloat, 3> &getColor() const { return _color; }
            std::array<GLfloat, 3> getLightColor() const
            {
                std::array<GLfloat, 3> retval = {{_color[0] * _intensity,
                                                  _color[1] * _intensity,
                                                  _color[2] * _intensity}};
                return retval;
            }

            /*! \brief Returns a projected light position.
     */
            magnet::math::Vector getEyespacePosition(const magnet::GL::Camera &camera) const
            {
                magnet::math::Vector vec = getPosition();
                magnet::math::NVector<GLfloat, 4> lightPos = {{GLfloat(vec[0]), GLfloat(vec[1]), GLfloat(vec[2]), 1.0f}};
                magnet::math::NVector<GLfloat, 4> lightPos_eyespace = camera.getViewMatrix() * lightPos;
                return magnet::math::Vector{lightPos_eyespace[0], lightPos_eyespace[1], lightPos_eyespace[2]};
            }
        };
    }
}
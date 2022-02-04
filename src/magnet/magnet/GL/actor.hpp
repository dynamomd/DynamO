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

namespace magnet
{
    namespace GL
    {
        namespace detail
        {
            enum class RenderMode
            {
                DEFAULT = 1 << 0, //!< The object is to render the standard data
                SHADOW = 1 << 1,  //!< This is a shadow pass (for lighting calculations).
                PICKING = 1 << 2  //!< This is an object picking pass.
            };
        }
        
        class Light;

        class Actor
        {
        protected:
            bool _visible = true;
            bool _shadow_casting = false;

        public:
            Actor() {}
            virtual ~Actor() {}
            virtual void deferred_pass(const Camera &, const detail::RenderMode){};
            virtual void forward_pass(const Camera &, const std::vector<std::shared_ptr<Light> >&, GLfloat,  const detail::RenderMode) {};

            bool visible() const { return _visible; }
            bool shadowCasting() const { return _shadow_casting; }
        };
    }
}
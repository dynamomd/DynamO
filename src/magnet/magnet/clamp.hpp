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
#include <algorithm>

namespace magnet {
  /*! \brief A C++ equivalent of the OpenCL and OpenGL clamp function.
   *
   * \param Val The value to clamp.
   * \param Min The lower bound to return if Val < Min.
   * \param Max The upper bound to return if Val > Max.
   */
  template<class T>
  inline T clamp(T Val, T Min, T Max) { return std::max(std::min(Val, Max), Min); }
}

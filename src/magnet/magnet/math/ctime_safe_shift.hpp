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

namespace magnet {
  namespace math {
    /*! \brief A template metafunction to perform a left bitwise shift.
     
      This function seems unnecessary at first, until you realize
      that a bit shift which is larger than the shifted type is an
      undefined operation. This function yields a sensible 0, instead
      of performing a partial shift, allowing template metaprograms
      to perform out of range shifts correctly.
     */    
    template<class T, T val, size_t shift, class Enable = void>
    struct ctime_safe_lshift {
      static const T result = val << shift;
    };

    template<class T, T val, size_t shift>
    struct ctime_safe_lshift<T, val, shift, typename std::enable_if<shift >= (sizeof(T) * 8)>::type> {
      static const T result = 0;
    };
    
    /*! \brief A template metafunction to perform a right bitwise shift.
     *
     * Please see \ref ctime_safe_lshift for the rational behind this
     * function.
     *
     * \sa ctime_safe_lshift
     */    
    template<class T, T val, size_t shift, class Enable = void>
    struct ctime_safe_rshift {
      static const T result = val >> shift;
    };

    template<class T, T val, size_t shift>
    struct ctime_safe_rshift<T, val, shift, typename std::enable_if<shift >= (sizeof(T) * 8)>::type> {
      static const T result = 0;
    };
  }
}


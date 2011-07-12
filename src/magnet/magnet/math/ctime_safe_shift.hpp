/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
//GCC fasely detects out of bounds shifts, even though they're
//protected behind the trinary condition! This just disables warnings
//for this small segment of code.
#pragma GCC system_header

namespace magnet {
  namespace math {
    /*! \brief A template metafunction to perform a left bitwise shift.
     *
     * This function seems unnessacary at first, until you realize
     * that a bit shift which is larger than the shifted type is an
     * undefined operation. This function yeilds a sensible 0, instead
     * of performing a partial shift, allowing template metaprograms
     * to perform out of range shifts correctly.
     */    
    template<class T, T val, size_t shift>
    struct ctime_safe_lshift {
      static const T result = (shift >= (sizeof(T) * 8)) ?  0 : (val << shift);
    };

    /*! \brief A template metafunction to perform a right bitwise shift.
     *
     * Please see \ref ctime_safe_lshift for the rational behind this
     * function.
     *
     * \sa ctime_safe_lshift
     */    
    template<class T, T val, size_t shift>
    struct ctime_safe_rshift {
      static const T result = (shift >= (sizeof(T) * 8)) ?  0 : (val >> shift);
    };
  }
}


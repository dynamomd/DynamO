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

namespace magnet {
  namespace math {
    //!\brief A template metafunction for calculating the power of an integer.
    template<int X, int Y>
    struct ctime_pow {
      static const int result = X * ctime_pow<X, Y-1>::result;
    };
    
    template<int X>
    struct ctime_pow<X,1> {
      static const int result = X;
    };

    template<int X>
    struct ctime_pow<X,0> {
      static const int result = 1;
    };
  }
}


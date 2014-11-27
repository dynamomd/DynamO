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
#include <cmath>
#include <limits>
#include <type_traits>

namespace magnet {
  namespace math {
    /*! \name Math functions
      \{
    */
    
    /*! \brief Calculate a "precision" for subtraction between two
      float types.
      
      If two floats are being subtracted, the precision of the
      operation is related to how close the numbers are. If they are
      of the same magnitude the precision may be terrible (so called
      catastrophic cancellation). We can rank how precise the
      operation is by comparing the difference in their exponents. A
      larger difference is always better, so the absolute difference
      is returned by this function.

      This function also handles the special cases where one or more
      of the arguments.
    */
    template<class T>
    size_t subtraction_precision(const T f1, const T f2) {
      static_assert(std::is_floating_point<T>(), "Can only calculate the precision of addition between floating point types");
      
      //Catch the case where this is not actually a subtraction at all
      if ((f1 == 0) || (f2 == 0) || (std::signbit(f1) != std::signbit(f2)))
	return std::numeric_limits<size_t>::max();
      
      int exp1, exp2;
      std::frexp(f1, &exp1);
      std::frexp(f2, &exp2);
      return std::abs(exp1-exp2);
    }

    /*! \brief Calculate a "precision" for addition between two float
      types.
      
      See subtraction_precision.

      This function also handles the special cases where one or more
      of the arguments.
    */
    template<class T>
    size_t addition_precision(const T f1, const T f2) {
      return subtraction_precision(f1, -f2);
    }
    
    /*! \} */
  }
}

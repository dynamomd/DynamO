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
#include <cstddef>

namespace magnet {
  namespace math {
    namespace detail {
      template <size_t val, size_t base>
      struct ctime_floor_log_worker
      {
	static const size_t result = ctime_floor_log_worker<val / base, base>::result + 1;
      };
      
      template <size_t base>
      struct ctime_floor_log_worker<0, base>
      {
	static const size_t result = 0;
      };

      template <size_t val, size_t base, size_t remainder>
      struct ctime_ceil_log_tester
      {
	static const size_t result = 1;
      };

      template <size_t base>
      struct ctime_ceil_log_tester<1,base,0>
      {
	static const size_t result = 0;
      };

      template <size_t base>
      struct ctime_ceil_log_tester<0,base,0>
      {
	static const size_t result = 0;
      };

      template <size_t val, size_t base>
      struct ctime_ceil_log_tester<val, base,0>
      {
	static const size_t result = ctime_ceil_log_tester<val / base, base, val % base>::result;
      };
    }

    /*! \brief A template metafunction to calculate the log of an
     * size_t integer.
     *
     * This function actually returns \f${\textrm
     * ceil}\left(\log_base(val)\right)\f$.
     */
    template <size_t val, size_t base>
    struct ctime_floor_log
    {
      static const size_t result = detail::ctime_floor_log_worker<val / base, base>::result;
    };
    
    // A specialization to produce an error for the log of zero.
    template <size_t base> struct ctime_floor_log<0, base> {};

    // A specialization to produce an error for the base 1 log.
    template <size_t val> struct ctime_floor_log<val, 1> {};

    // A specialization to produce an error for the base 0 log.
    template <size_t val> struct ctime_floor_log<val, 0> {};
    
    /*! \brief A template metafunction to calculate the log of an
     * size_t integer.
     *
     * This function actually returns \f${\textrm
     * floor}\left(\log_base(val)\right)\f$.
     */
    template <size_t val, size_t base>
    struct ctime_ceil_log
    {
      static const size_t result = ctime_floor_log<val,base>::result 
	+ detail::ctime_ceil_log_tester<val / base, base, val % base>::result;
    };

    template <size_t base>
    struct ctime_ceil_log<1, base>
    {
      static const size_t result = 0;
    };

    // A specialization to produce an error for the log of zero.
    template <size_t base> struct ctime_ceil_log<0, base> {};

    // A specialization to produce an error for the base 1 log.
    template <size_t val> struct ctime_ceil_log<val, 1> {};

    // A specialization to produce an error for the base 0 log.
    template <size_t val> struct ctime_ceil_log<val, 0> {};
  }
}

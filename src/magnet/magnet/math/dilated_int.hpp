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
/*! \file dilatedint.hpp
 * \brief Contains DilatedInteger and DilatedVector
 */

#pragma once

#include <stdint.h>
#include <limits>
#include <magnet/math/ctime_log.hpp>
#include <magnet/math/ctime_pow.hpp>

namespace magnet {
  namespace math {
    namespace detail {
      /*! \brief A class which provides all the constants needed for
       * dilated integers.
       */
      template<uint32_t _d>
      class DilatedConstants
      {
	/*! \brief Number of bits in the uint32_t type.
	 */
	static const uint32_t uint_bits = sizeof(uint32_t) * 8;

	/*! \brief Minimum number of usable bits in the dilated integer (\f$s\$).
	 *
	 * This is a minimum, as there may be 1 extra bit available for
	 * some values of the dilation. E.g., with a 32bit uint32_t, and a
	 * dilation of 3, you may interleave two 11-bit values and one
	 * 10-bit value to make a 32bit 3D Morton number.
	 *
	 * However, here we assume it is the maximum value!
	 */
	static const uint32_t _s = uint_bits / _d;
      
	/*! \brief A mask for the digit bits in the undilated integer.
	 */
	static const uint32_t undilatedMask = (uint32_t(0) - uint32_t(1)) >> (uint_bits - _s);

	/*! \brief A compile time calculation of \f${\textrm floor}(\log_{d-1}(s))\f$.
	 */
	static const uint32_t _t = ctime_log<_s, _d - 1>::value;

	/*! \brief A compile time calculation of \f${\textrm floor}(\log_d(s))\f$.
	 */
	static const uint32_t log_d_s = ctime_log<_s, _d>::value;

	//! \brief The calculator for the DilatedInteger x function.
	template <uint32_t p, uint32_t q>
	struct xworker
	{
	  static const uint32_t value = xworker<p-1,q>::value +  (1 << (p * q));
	};

	//! \brief The calculator for the DilatedInteger x function.
	template <uint32_t q>
	struct xworker<0,q>
	{
	  static const uint32_t value = 1;
	};

	template <uint32_t p, uint32_t q>
	struct x
	{
	  static const uint32_t value = xworker<p, q>::value;
	};

	template <uint32_t i_plus_1>
	struct c
	{
	  static const uint32_t value = x<_d, (_d-1) * ctime_pow<_d, i_plus_1 - 1>::result>::value;
	};

	template <uint32_t i>
	struct b
	{
	  static const uint32_t value = x<_d, ctime_pow<_d-1, _t - i + 1>::result>::value;
	};
      }
    }


    template<uint32_t _d>
    class DilatedInteger;

    template<>
    class DilatedInteger<1>; // 1-dilated objects are not dilated integers!
    
    /*! \brief A generic class for dilated integers.
     *
     * This class uses compile time integer mathematics to calculate
     * the constants for the dilation and undilation of the integers.
     *
     * All math symbols are taken from the excellent paper "Converting
     * to and from Dilated Integers", by R. Ramen and D. S. Wise
     * (doi:10.1109/TC.2007.70814). The mathematics on dilated
     * integers is explained by "Fast Additions on Masked
     * Integers"(10.1145/1149982.1149987).
     */
    template<uint32_t _d>
    class DilatedInteger
    {
    public:
      // Constructors
      inline DilatedInteger(): value(0) {}

      inline DilatedInteger(const uint32_t val):
	value(val) {}
      
      inline const uint32_t& getDilatedVal() const { return value; }
      
    private:
      uint32_t value; // stored as normalized integer at mask’s 1 bits.
    };   
  }
}

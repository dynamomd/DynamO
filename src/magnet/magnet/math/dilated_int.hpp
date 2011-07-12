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
#include <magnet/math/ctime_log.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/math/ctime_safe_shift.hpp>

namespace magnet {
  namespace math {
    namespace detail {
      /*! \brief A class which provides all the implementation details
       * needed for dilated integers.
       *
       * This class also contains all of the template metaprograms for
       * calculating the constants needed in dilating and undilating
       * these integers.
       */
      template<size_t _d>
      struct DilatedConstants
      {
	/*! \brief Number of bits in the size_t type.
	 */
	static const size_t uint_bits = sizeof(size_t) * 8;

	/*! \brief Minimum number of usable bits in the dilated integer (\f$s\$).
	 *
	 * This is a minimum, as there may be 1 extra bit available for
	 * some values of the dilation. E.g., with a 32bit size_t, and a
	 * dilation of 3, you may interleave two 11-bit values and one
	 * 10-bit value to make a 32bit 3D Morton number.
	 *
	 * However, here we assume it is the maximum value!
	 */
	static const size_t _s = uint_bits / _d;
      
	/*! \brief A mask for the digit bits in the undilated integer.
	 */
	static const size_t undilatedMask = ctime_safe_rshift<size_t,
							       (size_t(0) - size_t(1)), 
							       (uint_bits - _s)>::result;

	/*! \brief The number of rounds in a undilation.
	 *
	 * A compile time calculation of \f${\textrm ceil}(\log_{d}(s))\f$
	 */
	static const size_t undilation_rounds = ctime_ceil_log<_s, _d>::result; 

	/*! \brief The number of rounds in a dilation.
	 *
	 * A compile time calculation of \f${\textrm ceil}(\log_{d-1}(s))\f$
	 */
	static const size_t dilation_rounds = ctime_ceil_log<_s, _d - 1>::result;

	template <size_t p, size_t q>
	struct xworker
	{
	  static const size_t result = xworker<p-1, q>::result + ctime_safe_lshift<size_t, 1, p * q>::result;
	};

	template <size_t q>
	struct xworker<0,q>
	{
	  static const size_t result = 1;
	};

	//! \brief The calculator for the \f$x_{p,q}\f$ constant.
	template <size_t p, size_t q>
	struct x
	{
	  static const size_t result = xworker<p - 1, q>::result;
	};

	//! \brief The calculator for the \f$c_{d,i}\f$ constant.
	template <size_t i_plus_1>
	struct c
	{
	  static const size_t result = x<_d, (_d - 1) * ctime_pow<_d, i_plus_1 - 1>::result>::result;
	};

	//! \brief The calculator for the \f$b_{d,i}\f$ constant.
	template <size_t i>
	struct b
	{
	  static const size_t result = x<_d, ctime_pow<_d - 1, dilation_rounds - i + 1>::result>::result;
	};

	//! \brief Produces a \ref result with the lowest \ref n bits set.
	template <size_t n>
	struct getnbits
	{	  
	  static const size_t result = (ctime_safe_lshift<size_t,1, n>::result - 1);
	};

	template <size_t i, size_t counter>
	struct zworker
	{
	  static const size_t _di = ctime_pow<_d, i>::result;
	  static const size_t result = zworker<i, counter - 1>::result
	    | ctime_safe_rshift<size_t, zworker<i, 0>::result, ((_di * _d) * counter)>::result;
	};

	template <size_t i>
	struct zworker<i, 0>
	{
	  static const size_t _di = ctime_pow<_d, i>::result;
	  static const size_t _bitcount = (_di < _s) ? _di : _s;
	  static const size_t result = 
	    ctime_safe_lshift<size_t,
			      getnbits<_bitcount>::result & getnbits<_s>::result,
			      _d * (_s - 1) + 1 - _bitcount>::result;
	};

	//! \brief Calculator for the \f$z_{d,i}\f$ constant.
	template <size_t i>
	struct z
	{
	  static const size_t result = zworker<i, _s / ctime_pow<_d, i>::result>::result;
	};


	template <size_t i, size_t counter>
	struct yworker
	{
	  static const size_t _bitcount = ctime_pow<_d -1, dilation_rounds - i>::result;
	  static const size_t _bitseperation = ctime_pow<_d-1, dilation_rounds -i + 1>::result + _bitcount;
	  static const size_t result = yworker<i, counter - 1>::result 
	    | ctime_safe_lshift<size_t, getnbits<_bitcount>::result, _bitseperation * counter>::result;
	};

	template <size_t i>
	struct yworker<i, 0>
	{
	  static const size_t _bitcount = ctime_pow<_d -1, dilation_rounds - i>::result;
	  static const size_t result = getnbits<_bitcount>::result;
	};

	//! \brief Calculator for the \f$z_{d,i}\f$ constant.
	template <size_t i>
	struct y
	{
	  static const size_t _bitcount = ctime_pow<_d - 1, dilation_rounds - i>::result;
	  static const size_t result = yworker<i, _s / _bitcount>::result;
	};


	static const size_t undilate_shift = (_d * (_s - 1) + 1 - _s); 
	template <size_t i, int FakeParam>
	struct undilateWorker
	{
	  static inline size_t eval(size_t val)
	  {
	    return (undilateWorker<i - 1, FakeParam>::eval(val) * c<i>::result) & z<i>::result;
	  }
	};

	template <int FakeParam>
	struct undilateWorker<1, FakeParam>
	{
	  static inline size_t eval(size_t val)
	  {
	    return (val * c<1>::result) & z<1>::result;
	  }
	};
	
	//! \brief The undilation function.
	static inline size_t undilate(size_t val)
	{
	  return undilateWorker<undilation_rounds, 0>::eval(val) >> undilate_shift;
	};

	template <size_t i, int FakeParam>
	struct dilateWorker
	{
	  static inline size_t eval(size_t val)
	  {
	    return (dilateWorker<i - 1, FakeParam>::eval(val) * b<i>::result) & y<i>::result;
	  }
	};

	template <int FakeParam>
	struct dilateWorker<1, FakeParam>
	{
	  static inline size_t eval(size_t val)
	  {
	    return (val * b<1>::result) & y<1>::result;
	  }
	};
	
	//! \brief The dilation function.
	static inline size_t dilate(size_t val)
	{
	  return dilateWorker<dilation_rounds, 0>::eval(val);
	};
      };
    }


    template<size_t _d>
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
    template<size_t _d>
    class DilatedInteger
    {
    public:
      // Constructors
      inline DilatedInteger(): value(0) {}

      inline DilatedInteger(const size_t val):
	value(val) {}
      
      inline const size_t& getDilatedVal() const { return value; }
      
    private:
      size_t value; // stored as normalized integer at mask’s 1 bits.
    };   
  }
}

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
#include <stdexcept>

namespace magnet {
  namespace math {
    namespace dilatedinteger {
      /*! \brief Number of bits in the size_t type.
       */
      static const size_t uint_bits = sizeof(size_t) * 8;
      
      /*! \brief The number of usable bits in the dilated integer
       * (\f$s\$).
       *
       * This is technically a minimum value, as there may be 1 extra
       * bit available for some values of the dilation. E.g., with a
       * 32bit size_t, and a dilation of 3, you may interleave two
       * 11-bit values and one 10-bit value to make a 32bit 3D Morton
       * number.
       *
       * But please, don't try to use that extra bit, just use a 64
       * bit computer instead.
       */
      template <size_t d>
      struct s
      { static const size_t result = uint_bits / d; };
      
      /*! \brief The number of rounds in a dilation.
       *
       * A compile time calculation of \f${\textrm ceil}(\log_{d-1}(s))\f$
       */
      template <size_t d>
      struct dilation_rounds
      { static const size_t result  = ctime_ceil_log<s<d>::result, d - 1>::result; };

      template <>
      struct dilation_rounds<2>
      { static const size_t result  = ctime_ceil_log<s<2>::result, 2>::result; };

      /*! \brief The number of rounds in a undilation.
       *
       * A compile time calculation of \f${\textrm ceil}(\log_{d}(s))\f$
       */
      template <size_t d>
      struct undilation_rounds
      { static const size_t result  = ctime_ceil_log<s<d>::result, d>::result; };

      //! \brief The calculator for the \f$x_{p,q}\f$ constant.
      template <size_t p, size_t q>
      struct x
      {
	template <size_t l, size_t dummy>
	struct xworker
	{ static const size_t result = xworker<l - 1, dummy>::result + ctime_safe_lshift<size_t, 1, l * q>::result; };

	template <size_t dummy>
	struct xworker<0, dummy>
	{ static const size_t result = 1; };
	
	static const size_t result = xworker<p - 1, 0>::result;
      };

      //! \brief The calculator for the \f$c_{d,i}\f$ constant.
      template <size_t i_plus_1, size_t d>
      struct c
      { static const size_t result = x<d, (d - 1) * ctime_pow<d, i_plus_1 - 1>::result>::result; };

      //! \brief The calculator for the \f$b_{d,i}\f$ constant.
      template <size_t i, size_t d>
      struct b
      { static const size_t result = x<d, ctime_pow<d - 1, dilation_rounds<d>::result - i + 1>::result>::result; };
      
      //! \brief Produces a \ref result with the lowest \ref n bits set.
      template <size_t n>
      struct getnbits
      { static const size_t result = (ctime_safe_lshift<size_t,1, n>::result - 1); };
      
      //! \brief Returns the maximum value that can be dilated
      template <size_t d>
      struct maxDilatableValue
      { static const size_t result = getnbits<s<d>::result>::result; };

      //! \brief Calculator for the \f$z_{d,i}\f$ constant.
      template <size_t i, size_t d>
      struct z
      {
	static const size_t _di = ctime_pow<d, i>::result;
	static const size_t _bitcount = (_di < s<d>::result) ? _di : s<d>::result;

	template <size_t counter, size_t dummy>
	struct zworker
	{
	  static const size_t result = zworker<counter - 1, dummy>::result
	    | ctime_safe_rshift<size_t, zworker<0, dummy>::result, ((_di * d) * counter)>::result;
	};
	
	template <size_t dummy>
	struct zworker<0, dummy>
	{
	  static const size_t result = 
	    ctime_safe_lshift<size_t,
			      getnbits<_bitcount>::result & getnbits<s<d>::result>::result,
			      d * (s<d>::result - 1) + 1 - _bitcount>::result;
	};
	
	static const size_t result = zworker<s<d>::result / ctime_pow<d, i>::result, 0>::result;
      };
      
      //! \brief Calculator for the \f$z_{d,i}\f$ constant.
      template <size_t i, size_t d>
      struct y
      {
	static const size_t _bitcount = ctime_pow<d - 1, dilation_rounds<d>::result - i>::result;
	static const size_t _bitseperation = ctime_pow<d - 1, dilation_rounds<d>::result - i + 1>::result 
	  + _bitcount;

	template <size_t counter, size_t dummy>
	struct yworker
	{
	  static const size_t result = yworker<counter - 1, dummy>::result 
	    | ctime_safe_lshift<size_t, getnbits<_bitcount>::result, _bitseperation * counter>::result;
	};
	
	template <size_t dummy>
	struct yworker<0, dummy>
	{
	  static const size_t result = getnbits<_bitcount>::result;
	};
	
	static const size_t result = yworker<s<d>::result / _bitcount, 0>::result;
      };

      //! \brief The undilation function.
      template<size_t d> 
      struct undilate
      {
	template <size_t i, size_t dummy>
	  struct undilateWorker
	{
	  static inline size_t eval(size_t val)
	  {
	    return (undilateWorker<i - 1, dummy>::eval(val) * c<i, d>::result) & z<i, d>::result;
	  }
	};
	
	template <size_t dummy>
	struct undilateWorker<1, dummy>
	{
	  static inline size_t eval(size_t val)
	  {
	    return (val * c<1, d>::result) & z<1, d>::result;
	  }
	};

	static inline size_t eval(size_t val)
	{
	  static const size_t undilate_shift = (d * (s<d>::result - 1) + 1 - s<d>::result); 
	  return undilateWorker<undilation_rounds<d>::result, 0>::eval(val) >> undilate_shift;
	};
      };

      //! \brief The dilation function.
      template<size_t d> 
      struct dilate
      {
	template <size_t i, size_t dummy>
	struct dilateWorker
	{
	  static inline size_t eval(size_t val)
	  { return (dilateWorker<i - 1, d>::eval(val) *  b<i,d>::result) & y<i,d>::result; }
	};
	
	template <size_t dummy>
	struct dilateWorker<1, dummy>
	{
	  static inline size_t eval(size_t val)
	  { return (val * b<1, d>::result) & y<1, d>::result; }
	};

	static inline size_t eval(size_t val)
	{
	  return dilateWorker<dilation_rounds<d>::result, 0>::eval(val);
	};
      };

      /*! \brief A specialization of dilate for d=2.
       *
       * We have to use a specialization of dilate for a dilation
       * width of \f$d==2\f$, as the multiplication method is not
       * valid here.  Instead, we use the so-called Shift-Or
       * algorithm.
       */
      template<>
      struct dilate<2>
      {
	template <size_t i>
	struct y
	{
	  static const size_t _bitcount = s<2>::result / ( 1 << i);
	  static const size_t _bitseperation = 2 * _bitcount;
	  
	  template <size_t counter, size_t dummy>
	  struct yworker
	  {
	    static const size_t result = yworker<counter - 1, dummy>::result 
	      | ctime_safe_lshift<size_t, getnbits<_bitcount>::result, _bitseperation * counter>::result;
	  };
	  
	  template <size_t dummy>
	  struct yworker<0, dummy>
	  {
	    static const size_t result = getnbits<_bitcount>::result;
	  };
	  
	  static const size_t result = yworker<s<2>::result / _bitcount, 0>::result;
	};

	
	template <size_t i>
	struct shiftval
	{ static const size_t result = ctime_safe_lshift<size_t, 1, dilation_rounds<2>::result - i>::result; };

	template <size_t i, size_t dummy>
	struct dilateWorker
	{
	  static inline size_t eval(size_t val)
	  { 
	    const size_t val2 = dilateWorker<i - 1, 2>::eval(val);
	    return (val2 | (val2 << shiftval<i>::result)) & y<i>::result;
	  }
	};

	template <size_t dummy>
	struct dilateWorker<1, dummy>
	{
	  static inline size_t eval(size_t val)
	  {
	    return (val | (val << shiftval<1>::result)) & y<1>::result;
	  }
	};

	static inline size_t eval(size_t val)
	{
	  return dilateWorker<dilation_rounds<2>::result, 0>::eval(val);
	};
      };
    }
    
    /*! \brief A function to dilate an integer value.
     *
     * Dilation is the act of spreading the bits of an integer out
     * over a wider range by introducing \ref d zero bits in between
     * every bit of the original integer. This amount of spread is
     * known as the dilation.
     *
     * This function accepts and returns size_t values as it is likely
     * you will use the dilated integer for memory access (Dilated
     * integers are required for the calculation of Morton numbers, a
     * nice way to increase cache efficiency).
     *
     * A design choice was made to dilated the integers using
     * mathematical/bitwise operations instead of using look-up
     * tables. This is because calculation is relatively cheap in
     * modern CPU's. Also, if you are using dilated integers to
     * optimize your memory access patterns, then we don't want to
     * fill the cache up with even more data.
     *
     * Please see the original paper "Converting to and from Dilated
     * Integers"(doi:10.1109/TC.2007.70814) for the details on the
     * underlying math. The rest is template metaprogramming.
     *
     * \tparam d The requested dilation of the passed integer.
     */
    template<size_t d>
    inline size_t dilate(size_t val)
    {
#ifdef MAGNET_DYNAMO
      if (val > dilatedinteger::maxDilatableValue<d>::result)
	throw << std::out_of_range("Cannot dilate value");
#endif
      return dilatedinteger::dilate<d>::eval(val);
    };

    /*! \brief A function to dilate an integer value.
     *
     * \note Please see \ref dilate for more information on what dilation is!
     * \sa dilate
     * \tparam d The dilation of the passed integer.
     */
    template<size_t d>
    inline size_t undilate(size_t val)
    {
      return dilatedinteger::undilate<d>::eval(val);
    };
  }
}

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
#include <magnet/math/ctime_log.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/math/ctime_safe_shift.hpp>
#include <stdexcept>
#include <limits>

/*! \file dilated_int.hpp
 *
 * \brief An implementation of functions for arbitrary dilation of size_t
 * types.
 *
 * This file contains two functions for dilating and undilating size_t
 * integers. they are used like so
 * \code //Performing a 3-dilation of an integer
 * size_t val = 10;
 * dilatedval = magnet::math::dilate<3>(val);
 * undilatedval = magnet::math::undilate<3>(val);
 * assert(undilatedval == val);
 * \endcode
 *
 * The code for dilating the size_t type is generated at compile time,
 * depending on the dilation you use. The dilated/undilated values are
 * passed/returned as size_t types, as this is the type used for
 * addressing memory. 
 
 * This means that if you have a 32/64bit architecture, then you will
 * generate 32/64bit dilated integers.
 *
 * Finally, there is a helper class (\ref magnet::math::DilatedInteger), which
 * represents a DilatedInteger to help you make type-safe code. It
 * also supports some simple math operations directly on the dilated
 * integer.
 *
 * All of the implementation details of the dilation constants are at
 * the top of this file. This is implemented using template
 * meta-programming and is a little abstract. So I recommend that if
 * you want to learn how this code works, I suggest you read the
 * document backwards (bottom to top)! 
 */
namespace magnet {
  namespace math {
    //! \brief Implementation details  for the \ref dilate and \ref undilate functions.
    namespace dilatedinteger {
#ifndef DOXYGEN_SHOULD_IGNORE_THIS
      /*! \brief Number of bits in the size_t type.
       */
      static const size_t uint_bits = sizeof(size_t) * 8;
      
      /*! \brief The number of usable bits in the dilated integer
       * (\f$s\f$).
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
      
      /*! \brief Produces a result with the lowest n bits set.
       *
       * \tparam n The number of set bits to generate.
       */
      template <size_t n>
      struct getnbits
      { static const size_t result = (ctime_safe_lshift<size_t,1, n>::result - 1); };
      
      //! \brief Returns the maximum value that can be dilated.
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
      
      /*! \brief Calculator for the \f$y_{d,i}\f$ constant.
       *
       * The \f$y_{d,i}\f$ constant is the bit mask used after each
       * round of the dilation algorithm.
       */
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

      /*! \brief Calculator for the \f$y_{d,i}\f$ constant when \f$d==2\f$.
       */
      template <size_t i>
      struct y<i, 2>
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

      /*! \brief Returns the maximum value in its dilated form.
       *
       * This value is simply the bit mask from the final round of the
       * dilation algorithm.
       */
      template <size_t d>
      struct maxDilatedValue
      { static const size_t result = y<dilation_rounds<d>::result, d>::result; };


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

      //! \brief The actual dilation function.
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
       
        We have to use a specialization of dilate for a dilation
        width of \f$d==2\f$, as the multiplication method is not
        valid here.  Instead, we use the so-called Shift-Or
        algorithm.
       */
      template<>
      struct dilate<2>
      {
	template <size_t i>
	struct shiftval
	{ static const size_t result = ctime_safe_lshift<size_t, 1, dilation_rounds<2>::result - i>::result; };

	template <size_t i, size_t dummy>
	struct dilateWorker
	{
	  static inline size_t eval(size_t val)
	  { 
	    const size_t val2 = dilateWorker<i - 1, 2>::eval(val);
	    return (val2 | (val2 << shiftval<i>::result)) & y<i, 2>::result;
	  }
	};

	template <size_t dummy>
	struct dilateWorker<1, dummy>
	{
	  static inline size_t eval(size_t val)
	  {
	    return (val | (val << shiftval<1>::result)) & y<1, 2>::result;
	  }
	};

	static inline size_t eval(size_t val)
	{
	  return dilateWorker<dilation_rounds<2>::result, 0>::eval(val);
	};
      };
#endif
    }
    
    /*! \brief A function to dilate an integer value.
     
      Dilation is the act of spreading the bits of an integer out
      over a wider range by introducing \ref d zero bits in between
      every bit of the original integer. This amount of spread is
      known as the dilation.
     
      This function accepts and returns size_t values as it is likely
      you will use the dilated integer for memory access (Dilated
      integers are required for the calculation of Morton numbers, a
      nice way to increase cache efficiency).
     
      A design choice was made to dilated the integers using
      mathematical/bitwise operations instead of using look-up
      tables. This is because calculation is relatively cheap in
      modern CPU's. Also, if you are using dilated integers to
      optimize your memory access patterns, then we don't want to
      fill the cache up with even more data.
     
      Please see the original paper "Converting to and from Dilated
      Integers"(doi:10.1109/TC.2007.70814) for the details on the
      underlying math. The rest is template metaprogramming.
     
      \tparam d The requested dilation of the passed integer.
     */
    template<size_t d>
    inline size_t dilate(size_t val)
    {
#ifdef MAGNET_DEBUG
      if (val > dilatedinteger::maxDilatableValue<d>::result)
	throw std::out_of_range("Cannot dilate such a large value");
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
#ifdef MAGNET_DEBUG
      if (val > dilatedinteger::maxDilatedValue<d>::result)
	throw std::out_of_range("Cannot undilate such a large value");
#endif
      return dilatedinteger::undilate<d>::eval(val);
    };


    /*! \brief A helper class which allows mathematics to be performed
      directly on a dilated integer.
     
      This class is based on the paper "Fast additions on masked
      integers" by M. D. Adams and D. S. Wise
      (doi:10.1145/1149982.1149987).
     
     */
    template<size_t d>
    class DilatedInteger
    {
    public:      
      /*! \brief Default constructor.
       */
      inline DilatedInteger() {}

      /*! \brief Construct a DilatedInteger from a undilated integer.
       */
      inline DilatedInteger(const size_t val):
	_value(dilate<d>(val)) {}
      
      /*! \brief Copy constructor.
       */      
      inline DilatedInteger(const DilatedInteger<d>& o):
	_value(o._value) {}
      
      /*! \brief Returns the dilated integer. */      
      inline const size_t& getDilatedValue() const { return _value; }
      /*! \brief Returns the undilated integer. */      
      inline size_t getRealValue() const { return undilate<d>(_value); }
      /*! \brief Sets this DilatedInteger to the passed dilated integer. */
      inline void setDilatedValue(const size_t& i) { _value = i & dilatedMask; }
      
      /*! \brief Assignment operator for undilated values. */
      inline void operator=(const size_t& i) { _value = dilate<d>(i); }
      /*! \brief Assignment operator. */
      inline void operator=(const DilatedInteger& i) { _value = i._value; }
      
      //! \brief Subtraction operator.
      inline DilatedInteger operator-(const DilatedInteger& o) const
      { return DilatedInteger((_value - o._value) & dilatedMask, 0); }
      
      //! \brief Addition operator.
      inline DilatedInteger operator+(const DilatedInteger& o) const
      { return DilatedInteger((_value + (~dilatedMask) + o._value) & dilatedMask, 0); }
      
      //! \brief Increment operator.
      inline DilatedInteger& operator++() 
      { _value = (_value - dilatedMask) & dilatedMask; return *this; }
      
      //! \brief Decrement operator.
      inline DilatedInteger& operator--() 
      { _value = (_value - 1) & dilatedMask; return *this; }

      //! \brief Inplace subtraction operator.
      inline DilatedInteger& operator-=(const DilatedInteger& o)
      { 
	_value -= o._value; 
	_value &= dilatedMask; 
	return *this;
      }

      //! \brief Inplace addition operator.
      inline DilatedInteger& operator+=(const DilatedInteger& o)
      { 
	_value += (~dilatedMask) + o._value;
	_value &= dilatedMask;
	return *this;
      }
      
      //! \brief Comparison operator.
      inline bool operator==(const DilatedInteger& o) const
      { return _value == o._value; }
      
      //! \brief Inequality operator.
      inline bool operator!=(const DilatedInteger& o) const
      { return _value != o._value; }
      
      //! \brief Less-than operator.
      inline bool operator<(const DilatedInteger& o) const
      { return _value < o._value; }
      
      //! \brief Greater-than operator.
      inline bool operator>(const DilatedInteger& o) const
      { return _value > o._value; }
      
      //! \brief Less-than-or-equal operator.
      inline bool operator<=(const DilatedInteger& o) const
      { return _value <= o._value; }
      
      //! \brief Greater-than-or-equal operator.
      inline bool operator>=(const DilatedInteger& o) const
      { return _value >= o._value; }

      /*! \brief Modulus operator (expensive).
       
        A complication with the modulus operation is that the real
        value actually has less bits than the type. This means that
        large overflows are not wrapped correctly.
       
        To fix this we simply mask off the unused bits of the real
        value.
       */
      inline DilatedInteger operator%(const size_t& mod) const
      { return DilatedInteger(getRealValue() % mod); }

      //! \brief Inplace modulus operator (expensive).
      inline DilatedInteger& operator%=(size_t mod)
      {
	operator=(getRealValue() % mod);
	return *this;
      }
      
    private:
      template<class T> friend struct std::numeric_limits; 

      /*! \brief A mask for the settable bits of the dilated integer.
       */
      static const size_t dilatedMask = dilatedinteger::maxDilatedValue<d>::result;

      /*! \brief Hidden constructor constructing from pre-dilated
       * integer values.
       *
       * This constructor is used by the numeric limits class
       */
      inline DilatedInteger(const size_t val, void*):
	_value(val) {}

      /*! \brief The dilated integer, stored as normalized integer at
       * mask’s 1 bits.
       */
      size_t _value; 
    };
  }
}


namespace std {
  /*! \brief A specialization of numeric_limits for
   * \ref magnet::math::DilatedInteger classes.
   */
  template <size_t d>
  struct numeric_limits<magnet::math::DilatedInteger<d> > 
  {
  public:
    static const bool is_specialized = true;
    
    inline static magnet::math::DilatedInteger<d> min() throw() 
    { return magnet::math::DilatedInteger<d>(0,0); }
    
    inline static magnet::math::DilatedInteger<d> max() throw()
    { return magnet::math::DilatedInteger<d>(magnet::math::dilatedinteger::maxDilatedValue<d>::result,0); }

    static const int  digits = magnet::math::dilatedinteger::s<d>::result;
    static const int  digits10 = (size_t(1) << digits) / 10;
    static const bool is_signed = false;
    static const bool is_integer = true;
    static const bool is_exact = true;
    static const int radix = 2;
    inline static magnet::math::DilatedInteger<d> epsilon() throw() 
    { return magnet::math::DilatedInteger<d>(1, 0); }
    
    static magnet::math::DilatedInteger<d> round_error() throw();
    
    static const int  min_exponent = 0;
    static const int  min_exponent10 = 0;
    static const int  max_exponent = 0;
    static const int  max_exponent10 = 0;
    
    static const bool has_infinity = false;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN = false;
    
    //static const float_denorm_style has_denorm = denorm absent;

    static const bool has_denorm_loss = false;
    
    static magnet::math::DilatedInteger<d> infinity() throw();
    static magnet::math::DilatedInteger<d> quiet_NaN() throw();
    static magnet::math::DilatedInteger<d> signaling_NaN() throw();
    static magnet::math::DilatedInteger<d> denorm_min() throw();
    
    static const bool is_iec559 = false;
    static const bool is_bounded = true;
    static const bool is_modulo = true;
    
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
  };
}

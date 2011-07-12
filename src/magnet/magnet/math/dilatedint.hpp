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
#include <magnet/math/dilated_int.hpp>

namespace magnet {
  namespace math {

    /*! \brief A class to allow performing math on a dilated integer.
     */
    template<size_t d>
    class DilatedInteger
    {
    public:
      //A mask for the dilated integer (also dilated max value)
      static const size_t dilatedMask = dilatedinteger::maxDilatedValue<d>::result;
      static const size_t binarydigits = dilatedinteger::s<d>::result;
      
      // Constructors
      inline DilatedInteger() {}

      inline DilatedInteger(const size_t val):
	value(dilate<d>(val)) {}
      
      inline DilatedInteger(const DilatedInteger<d>& o):
	value(o.value) {}
      
      inline const size_t& getDilatedVal() const { return value; }
      
      inline size_t getRealVal() const { return undilate<d>(value); }
      
      inline void setDilatedVal(const size_t& i) { value = i & dilatedMask; }
      
      inline void operator=(const size_t& i) { value = dilate<d>(i); }
      inline void operator=(const DilatedInteger& i) { value = i.value; }
      
      inline void zero() { value = 0; }
      
      inline bool isZero() const { return value == 0; }
      
      // Simple operators.
      inline DilatedInteger operator-(const DilatedInteger& o) const
      { return DilatedInteger((value - o.value) & dilatedMask, 0); }
      
      inline DilatedInteger operator+(const DilatedInteger& o) const
      { return DilatedInteger((value + (~dilatedMask) + o.value) & dilatedMask, 0); }
      
      inline DilatedInteger& operator++() 
      { value = (value - dilatedMask) & dilatedMask; return *this; }
      
      inline DilatedInteger& operator--() 
      { value = (value - 1) & dilatedMask; return *this; }

      inline DilatedInteger& operator-=(const DilatedInteger& o)
      { 
	value -= o.value; 
	value &= dilatedMask; 
	return *this;
      }

      inline DilatedInteger& operator+=(const DilatedInteger& o)
      { 
	value += (~dilatedMask) + o.value;
	value &= dilatedMask;
	return *this;
      }
      
      inline bool operator==(const DilatedInteger& o) const
      { return value == o.value; }
      
      inline bool operator!=(const DilatedInteger& o) const
      { return value != o.value; }
      
      inline bool operator<(const DilatedInteger& o) const
      { return value < o.value; }
      
      inline bool operator>(const DilatedInteger& o) const
      { return value > o.value; }
      
      inline bool operator<=(const DilatedInteger& o) const
      { return value <= o.value; }
      
      inline bool operator>=(const DilatedInteger& o) const
      { return value >= o.value; }

      
    private:
      template<class T> friend struct std::numeric_limits; 

      inline DilatedInteger(const size_t val, void*):
	value(val) {}

      size_t value; // stored as normalized integer at mask’s  1 bits.
    };
    
    template<size_t d>
    struct DilatedVector
    {
      inline DilatedVector() {}
      
      inline DilatedVector(const size_t& MortonNum)
      {
	for (size_t i(0); i < d; ++i)
	  data[i].setDilatedVal(MortonNum >> i);
      }
      
      inline DilatedVector(const size_t& x, const size_t& y, const size_t& z)
      {
	data[0] = x;
	data[1] = y;
	data[2] = z;
      }

      inline DilatedVector(const DilatedInteger<d>& x, 
			   const DilatedInteger<d>& y, 
			   const DilatedInteger<d>& z)
      {
	data[0] = x;
	data[1] = y;
	data[2] = z;
      }
      
      inline size_t getMortonNum()
      { 
	size_t retval = data[0].getDilatedVal();
	for (size_t i(1); i < d; ++i)
	  retval += data[i].getDilatedVal() << i;

	return retval;
      }
      
      inline DilatedVector operator+(const DilatedVector& o) const
      {
	DilatedVector retval;
	for (size_t i(0); i < d; ++i)
	  retval.data[i] = o.data[i] + data[i];
	
	return retval;
      }

      DilatedInteger<d> data[d];
    };
  }
}

namespace std {
  template <size_t d>
  struct numeric_limits<magnet::math::DilatedInteger<d> > 
  {
  public:
    static const bool is_specialized = true;
    
    inline static magnet::math::DilatedInteger<d> min() throw() 
    { return magnet::math::DilatedInteger<d>(0,0); }
    
    inline static magnet::math::DilatedInteger<d> max() throw() 
    { return magnet::math::DilatedInteger<d>(magnet::math::DilatedInteger<d>::dilatedMask,0); }

    static const int  digits = magnet::math::DilatedInteger<d>::binarydigits;
    static const int  digits10 = (1 << magnet::math::DilatedInteger<d>::binarydigits) / 10;
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

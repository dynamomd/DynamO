/*  DYNAMO:- Event driven molecular dynamics simulator 
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

namespace magnet {
  namespace math {
    
    /*! \brief A dilated integer class for easy morton order manipulations in 3D
     *
     * Please see the original papers "Converting to and from Dilated
     * Integers"(10.1109/TC.2007.70814) and "Fast Additions on Masked
     * Integers"(10.1145/1149982.1149987) for more information and a 2D implementation
     */
    class DilatedInteger
    {
    public:
      //Number of bits in the dilated integer (10 is max for 3 dilated ints)
      static const uint32_t digits = 10;
      
      //A mask for the number of bits in the dilated integer (also max value)
      static const uint32_t undilatedMask = 0xFFFFFFFF >> (32 - digits);
      
      //A mask for the dilated integer (also dilated max value)
      static const uint32_t dilatedMask = 0x49249249 & (0xFFFFFFFF >> (32 - 3 * digits));
      
      // Constructors
      inline DilatedInteger() {}

      inline DilatedInteger(const uint32_t val):
	value(dilate_3(val & undilatedMask)) {}
      
      //Constructor that takes the actual Dilated int as the arg
      inline DilatedInteger(const uint32_t val, void*):
	value(val) {}

      inline DilatedInteger(const DilatedInteger& d):
	value(d.value) {}
      
      inline const uint32_t& getDilatedVal() const { return value; }
      
      inline uint32_t getRealVal() const { return undilate_3(value); }
      
      inline void setDilatedVal(const uint32_t& i) { value = i & dilatedMask; }
      
      inline void operator=(const uint32_t& i) { value = dilate_3(i & undilatedMask); }
      inline void operator=(const DilatedInteger& i) { value = i.value; }
      
      inline void zero() { value = 0; }
      
      inline bool isZero() const { return value == 0; }
      
      // Simple operators.
      inline DilatedInteger operator-(const DilatedInteger& d) const
      { return DilatedInteger((value - d.value) & dilatedMask, 0); }
      
      inline DilatedInteger operator+(const DilatedInteger& d) const
      { return DilatedInteger((value + (~dilatedMask) + d.value) & dilatedMask, 0); }
      
      inline DilatedInteger& operator++() 
      { value = (value - dilatedMask) & dilatedMask; return *this; }
      
      inline DilatedInteger& operator--() 
      { value = (value - 1) & dilatedMask; return *this; }

      inline DilatedInteger& operator-=(const DilatedInteger& d)
      { 
	value -= d.value; 
	value &= dilatedMask; 
	return *this;
      }

      inline DilatedInteger& operator+=(const DilatedInteger& d)
      { 
	value += (~dilatedMask) + d.value;
	value &= dilatedMask;
	return *this;
      }
      
      inline bool operator==(const DilatedInteger& d) const
      { return value == d.value; }
      
      inline bool operator!=(const DilatedInteger& d) const
      { return value != d.value; }
      
      inline bool operator<(const DilatedInteger& d) const
      { return value < d.value; }
      
      inline bool operator>(const DilatedInteger& d) const
      { return value > d.value; }
      
      inline bool operator<=(const DilatedInteger& d) const
      { return value <= d.value; }
      
      inline bool operator>=(const DilatedInteger& d) const
      { return value >= d.value; }
      
    private:
      uint32_t value; // stored as normalized integer at mask’s  1 bits.
      
      //Could technically return a uint16_t
      inline uint32_t undilate_3(uint32_t t) const 
      {
	t = (t * 0x00015) & 0x0E070381;
	t = (t * 0x01041) & 0x0FF80001;
	t = (t * 0x40001) & 0x0FFC0000;
	return t >> 18;
      }
      
      //Could technically accept a uint16_t
      inline uint32_t dilate_3(uint32_t r) const
      {
	r =  (r * 0x10001) & 0xFF0000FF; 
	r =  (r * 0x00101) & 0x0F00F00F;
	r =  (r * 0x00011) & 0xC30C30C3;
	r =  (r * 0x00005) & 0x49249249;
	return r;
      }
      
    };
    
    struct DilatedVector
    {
      inline DilatedVector() {}
      
      inline DilatedVector(const uint32_t& MortonNum)
      {
	for (uint32_t i(0); i < 3; ++i)
	  data[i].setDilatedVal(MortonNum >> i);
      }
      
      inline DilatedVector(const uint32_t& x, const uint32_t& y, const uint32_t& z)
      {
	data[0] = x;
	data[1] = y;
	data[2] = z;
      }

      inline DilatedVector(const DilatedInteger& x, 
			   const DilatedInteger& y, 
			   const DilatedInteger& z)
      {
	data[0] = x;
	data[1] = y;
	data[2] = z;
      }
      
      inline uint32_t getMortonNum()
      { 
	return data[0].getDilatedVal() 
	  + (data[1].getDilatedVal() << 1) 
	  + (data[2].getDilatedVal() << 2); 
      }
      
      inline DilatedVector operator+(const DilatedVector& d) const
      {
	return DilatedVector(d.data[0] + data[0],
			     d.data[1] + data[1],
			     d.data[2] + data[2]);
      }

      DilatedInteger data[3];
    };
  }
}

namespace std {
  template <>
  struct numeric_limits<magnet::math::DilatedInteger> 
  {
  public:
    static const bool is_specialized = true;
    
    inline static magnet::math::DilatedInteger min() throw() 
    { return magnet::math::DilatedInteger(0,0); }
    
    inline static magnet::math::DilatedInteger max() throw() 
    { return magnet::math::DilatedInteger(magnet::math::DilatedInteger::dilatedMask,0); }

    static const int  digits = magnet::math::DilatedInteger::digits;
    static const int  digits10 = (1 << magnet::math::DilatedInteger::digits) / 10;
    static const bool is_signed = false;
    static const bool is_integer = true;
    static const bool is_exact = true;
    static const int radix = 2;
    inline static magnet::math::DilatedInteger epsilon() throw() 
    { return magnet::math::DilatedInteger(1,0); }
    
    static magnet::math::DilatedInteger round_error() throw();
    
    static const int  min_exponent = 0;
    static const int  min_exponent10 = 0;
    static const int  max_exponent = 0;
    static const int  max_exponent10 = 0;
    
    static const bool has_infinity = false;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN = false;
    
    //static const float_denorm_style has_denorm = denorm absent;

    static const bool has_denorm_loss = false;
    
    static magnet::math::DilatedInteger infinity() throw();
    static magnet::math::DilatedInteger quiet_NaN() throw();
    static magnet::math::DilatedInteger signaling_NaN() throw();
    static magnet::math::DilatedInteger denorm_min() throw();
    
    static const bool is_iec559 = false;
    static const bool is_bounded = true;
    static const bool is_modulo = true;
    
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
  };
}

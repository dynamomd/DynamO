/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

namespace magnet {
  namespace math {
    class DilatedInteger
    {
    public:
      static const uint32_t S = 10;
      static const uint32_t Smask = (((uint32_t)0) - 1) >> (sizeof(uint32_t)*CHAR_BIT - S);
      //Zeros are needed in the last two bits
      static const uint32_t mask  = 0x09249249;
      static const uint32_t maxVal = (((uint32_t)0) - 1) >> (sizeof(uint32_t)*CHAR_BIT - S);
      static const uint32_t DilatedMaxVal = ((((uint32_t)0) - 1) >> (sizeof(uint32_t)*CHAR_BIT - S * 3)) & mask; 
      
      // Constructor and Getter
      DilatedInteger() {}
      DilatedInteger(const uint32_t& val):
	value(dilate_3(val & Smask)) {}
      
      //Constructor that takes the actual Dilated int as the arg
      DilatedInteger(const uint32_t& val, void*):
	value(val) {}
      
      const uint32_t& getDilatedVal() const { return value; }
      
      uint32_t getRealVal() const { return undilate_3(value); }
      
      void setDilatedVal(const uint32_t& i) { value = i & mask; }
      
      void operator=(const uint32_t& i) { value = dilate_3(i & Smask); }
      
      void zero() { value = 0; }
      
      bool isZero() const { return value == 0; }
      
      // Simple operators.
      DilatedInteger operator-(const DilatedInteger& d) const
      { return DilatedInteger((value - d.value) & mask, 0); }
      
      DilatedInteger operator+(const DilatedInteger& d) const
      { return DilatedInteger((value + (~mask) + d.value) & mask, 0); }
      
      DilatedInteger& operator++() { value = (value - mask) & mask; return *this; }
      
      DilatedInteger& operator--() { value = (value -    1) & mask; return *this; }
      
      bool operator==(const DilatedInteger& d) const
      { return value == d.value; }
      
      bool operator!=(const DilatedInteger& d) const
      { return value != d.value; }
      
      bool operator<(const DilatedInteger& d) const
      { return value < d.value; }
      
      bool operator>(const DilatedInteger& d) const
      { return value > d.value; }
      
      bool operator<=(const DilatedInteger& d) const
      { return value <= d.value; }
      
      bool operator>=(const DilatedInteger& d) const
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
      DilatedVector() {}
      
      DilatedVector(const uint32_t& MortonNum)
      {
	for (uint32_t i(0); i < 3; ++i)
	  data[i].setDilatedVal(MortonNum >> i);
      }
      
      DilatedVector(const uint32_t& x, const uint32_t& y, const uint32_t& z)
      {
	data[0] = x;
	data[1] = y;
	data[2] = z;
      }
      
      inline uint32_t getMortonNum()
      { return data[0].getDilatedVal() + (data[1].getDilatedVal() << 1) + (data[2].getDilatedVal() << 2); }
      
      DilatedInteger data[3];
    };
    
//#include <iostream>
//    
//    template <typename T>
//    void print_bits ( T val, std::ostream& out )
//    {
//      T n_bits = sizeof ( val ) * CHAR_BIT;
//      
//      for ( unsigned i = 0; i < n_bits; ++i ) 
//    	{
//    	  out<< !!( val & 1 );
//    	  val >>= 1;
//    	}
//    }
  }
}

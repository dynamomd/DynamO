/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
 * \brief Contains a class to help with 3D Morton ordering
 */

#ifndef dilatedint_H
#define dilatedint_H
#include <limits.h>

class MaskedInteger 
{
public:
  static const unsigned int S = 10;
  static const unsigned int Smask;
  //Zeros are needed in the last two bits (hence the 0)
  static const unsigned int mask  = 0x09249249;
  static const unsigned int maxVal;
  static const unsigned int dilatedMaxVal; 

  // Constructor and Getter
  MaskedInteger() {}
  MaskedInteger(const unsigned int& val):
    value(dilate_3(val & Smask)) {}

  //Constructor that takes the actual dilated int as the arg
  MaskedInteger(const unsigned int& val, void*):
    value(val) {}
  
  const unsigned int& getDilatedVal() const { return value; }
  
  unsigned int getRealVal() { return undilate_3(value); }
  
  void setDilatedVal(const unsigned int& i) { value = i & mask; }

  void operator=(const unsigned int& i) { value = dilate_3(i & Smask); }

  void zero() { value = 0; }
  
  bool isZero() const { return value == 0; }

  // Simple operators.
  MaskedInteger operator-(const MaskedInteger& d) const
  { return MaskedInteger((value - d.value) & mask, 0); }

  MaskedInteger operator+(const MaskedInteger& d) const
  { return MaskedInteger((value + (~mask) + d.value) & mask, 0); }

  MaskedInteger& operator++() { value = (value - mask) & mask; return *this; }

  MaskedInteger& operator--() { value = (value -    1) & mask; return *this; }

  bool operator==(const MaskedInteger& d) const
  { return value == d.value; }

  bool operator!=(const MaskedInteger& d) const
  { return value != d.value; }

  bool operator<(const MaskedInteger& d) const
  { return value < d.value; }

  bool operator>(const MaskedInteger& d) const
  { return value > d.value; }

  bool operator<=(const MaskedInteger& d) const
  { return value <= d.value; }

  bool operator>=(const MaskedInteger& d) const
  { return value >= d.value; }

private:
  unsigned int value; // stored as normalized integer at mask’s  1 bits.

  inline unsigned int undilate_3(unsigned int t) const 
  {
    t = (t * 0x00015) & 0x0E070381;
    t = (t * 0x01041) & 0x0FF80001;
    t = (t * 0x40001) & 0x0FFC0000;
    return t >> 18;
  }

  inline unsigned int dilate_3(unsigned int r) const
  {
    r =  (r * 0x10001) & 0xFF0000FF; 
    r =  (r * 0x00101) & 0x0F00F00F;
    r =  (r * 0x00011) & 0xC30C30C3;
    r =  (r * 0x00005) & 0x49249249;
    return r;
  }

};

struct dilatedCoords
{
  dilatedCoords() {}

  dilatedCoords(const unsigned int& MortonNum)
  {
    for (unsigned int i(0); i < 3; ++i)
      data[i].setDilatedVal(MortonNum >> i);
  }

  dilatedCoords(const unsigned int& x, const unsigned int& y, const unsigned int& z)
  {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  inline unsigned int getMortonNum()
  { return data[0].getDilatedVal() + (data[1].getDilatedVal() << 1) + (data[2].getDilatedVal() << 2); }

  MaskedInteger data[3];
};

typedef MaskedInteger MI;

#include <iostream>

template <typename T>
void print_bits ( T val, std::ostream& out )
{
  T n_bits = sizeof ( val ) * CHAR_BIT;

  for ( unsigned i = 0; i < n_bits; ++i ) 
    {
      out<< !!( val & 1 );
      val >>= 1;
    }
}

#endif

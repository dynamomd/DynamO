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
/*! \file threadpool.hpp
 * \brief Contains the definition of CThreadPool
 */

#pragma once
#ifndef dilatedint_H
#define dilatedint_H
#include <limits.h>

template<typename T,
	 T mask, 
	 T dimension,
	 T S = 10>
class MaskedInteger{
  static const T Smask = (T(0) - 1) >> (sizeof(T)*CHAR_BIT - S);
public:
  // Constructor and Getter
  MaskedInteger() {}
  MaskedInteger(const T& val):
    value(dilate_3(val & Smask) << dimension) {}

  //Constructor that takes the actual dilated int as the arg
  MaskedInteger(const T& val, void*):
    value(val) {}
  
  const T& getDilatedVal() { return value; }
  
  T getRealVal() { return undilate_3(value >> dimension); }
  
  void setDilatedVal(const T& i) { value = i & mask; }

  void operator=(const T& i) { value = dilate_3(i & Smask) << dimension; }
  
  // Simple operators.
  MaskedInteger operator-(const MaskedInteger& d) const
  { return MaskedInteger((value - d.value) & mask, 0); }

  MaskedInteger operator+(const MaskedInteger& d) const
  { return MaskedInteger((value + (~mask) + d.value) & mask, 0); }

  void operator++() { value = (value - mask) & mask; }

  void operator--() { value = (value -    1) & mask; }

  void operator==(const MaskedInteger& d) const
  { return value == d.value; }

  void operator!=(const MaskedInteger& d) const
  { return value != d.value; }

  void operator<(const MaskedInteger& d) const
  { return value < d.value; }

  void operator>(const MaskedInteger& d) const
  { return value > d.value; }

  void operator<=(const MaskedInteger& d) const
  { return value <= d.value; }

  void operator>=(const MaskedInteger& d) const
  { return value >= d.value; }

private:
  T value; // stored as normalized integer at mask’s  1 bits.

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

typedef MaskedInteger<unsigned int, 0x49249249, 0> MI0;
typedef MaskedInteger<unsigned int, 0x92492492, 1> MI1;
typedef MaskedInteger<unsigned int, 0x24924924, 2> MI2;

struct dilatedCoords
{
  dilatedCoords(const unsigned int& MortonNum)
  {
    _x.setDilatedVal(MortonNum);
    _y.setDilatedVal(MortonNum);
    _z.setDilatedVal(MortonNum);
  }

  dilatedCoords(const unsigned int& x, const unsigned int& y, const unsigned int& z):
    _x(x), _y(y), _z(z)
  {}

  inline unsigned int getMortonNum()
  { return _x.getDilatedVal() + _y.getDilatedVal() + _z.getDilatedVal(); }

  MI0 _x;
  MI1 _y;
  MI2 _z;
};

/*
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

int main()
{
  dilatedCoords test(1023,3,0);
  
  std::cout << test.getMortonNum() << std::endl;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  ++test._x;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  --test._y;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  test._y = test._y + MI1(5);
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  --test._z;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  test = dilatedCoords(153391707);
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  test._x = 21;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;
  
  std::cout << std::endl;

  print_bits(MI0(1).getDilatedVal(),std::cout);
  std::cout << std::endl;

  print_bits((MI0(2)).getDilatedVal(),std::cout);
  std::cout << std::endl;

  MI0 i(1);
  print_bits((i+2).getDilatedVal(), std::cout);
  std::cout << std::endl;

  std::cout << (i+2).getRealVal() << std::endl;
}
*/
#endif

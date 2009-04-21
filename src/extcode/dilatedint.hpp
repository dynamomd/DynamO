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

template<typename T,
	 T mask, 
	 T dimension,
	 T S = 10,
	 T Smask = 0x03FF>
class MaskedInteger{
public:
    // Constructor and Getter
  MaskedInteger(const T& val): 
    value(dilate_3(val & Smask) << dimension) {}

  inline const T& getDilatedVal() { return value; }

  inline T getRealVal() { return undilate_3(value >> dimension); }

  // Simple operators.
  MaskedInteger operator-(const MaskedInteger& d) const
  { return MaskedInteger((value - d.value) & mask); }
  MaskedInteger operator+(const MaskedInteger& d) const
  { return MaskedInteger((value + (~mask) + d.value) & mask); }
  inline void operator++() { value = (value - mask) & mask; }
  inline void operator--() { value = (value -    1) & mask; }
  inline void operator==(const MaskedInteger& d) const
  { return value == d.value; }
  inline void operator!=(const MaskedInteger& d) const
  { return value != d.value; }
  inline void operator<(const MaskedInteger& d) const
  { return value < d.value; }
  inline void operator>(const MaskedInteger& d) const
  { return value > d.value; }
  inline void operator<=(const MaskedInteger& d) const
  { return value <= d.value; }
  inline void operator>=(const MaskedInteger& d) const
  { return value >= d.value; }

private:
  T value; // stored as normalized integer at mask’s  1 bits.

  inline unsigned short undilate_3(unsigned int t){
    t = (t * 0x00015) & 0x0E070381;
    t = (t * 0x01041) & 0x0FF80001;
    t = (t * 0x40001) & 0x0FFC0000;
    return ((unsigned short) (t >> 18));
  }

  inline unsigned int dilate_3(unsigned short t) {
    unsigned int r = t;
    r =  (r * 0x10001) & 0xFF0000FF; 
    r =  (r * 0x00101) & 0x0F00F00F;
    r =  (r * 0x00011) & 0xC30C30C3;
    r =  (r * 0x00005) & 0x49249249;
    return(r);
  }

};

typedef MaskedInteger<unsigned int, 0x49249249, 0> MI0;
typedef MaskedInteger<unsigned int, 0x92492492, 1> MI1;
typedef MaskedInteger<unsigned int, 0x24924924, 2> MI2;

struct dilatedCoords 
{
  dilatedCoords(const unsigned short& x, const unsigned short& y, const unsigned short& z):
    _x(x), _y(y), _z(z)
  {}

  inline unsigned int getMortonNum()
  { return _x.getDilatedVal() + _y.getDilatedVal() + _z.getDilatedVal(); }

  MI0 _x;
  MI1 _y;
  MI2 _z;
};

/*int main()
{
  dilatedCoords test(1023,1023,1023);
  
  ++test._x;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  ++test._y;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;

  ++test._z;
  std::cout << test._x.getRealVal() << " " << test._y.getRealVal() << " " << test._z.getRealVal() <<  std::endl;
  
}
*/
#endif

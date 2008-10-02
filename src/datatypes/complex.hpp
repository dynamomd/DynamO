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

#ifndef CComplex_H
#define CComplex_H

#include <math.h>

class CComplex
{
 public:
  CComplex():r(0.0), i(0.0) {};
  CComplex(Iflt r2, Iflt i2):r(r2), i(i2) {};
  
  CComplex &operator+=(const CComplex &cc)
    {
      r += cc.r;
      i += cc.i;
      return *this;
    }

  CComplex &operator-=(const CComplex &cc)
    {
      r -= cc.r;
      i -= cc.i;
      return *this;
    }

  CComplex operator*(const CComplex &cc)
    {
      return CComplex(r*cc.r-i*cc.i, r*cc.i + i*cc.r);
    }

  CComplex operator/(const CComplex &cc)
    {
      Iflt denom = cc.r * cc.r + cc.i * cc.i;
      return CComplex((r * cc.r + i * cc.i)/denom, 
		      (i * cc.r - r * cc.i)/denom);
    }

  CComplex operator*(const Iflt &a)
    {
      return CComplex(r*a, i*a);
    }

  CComplex operator/(const Iflt &a)
    {
      return CComplex(r/a, i/a);
    }


  bool operator==(const CComplex &cc) const
    {
      return ((cc.r == r) && (cc.i == i));
    }

  Iflt modulus() const
    {
      Iflt tmp;
      if ((tmp = r*r + i*i) != 0.0)
	return sqrt(tmp);
      else
	return 0.0;
    }

  CComplex conjugate()
    {
      return CComplex(r,-i);
    }

  CComplex exponent()
    {
      return CComplex(exp(r)*cos(i),exp(r)*sin(i));
    }

  const Iflt &geti() {return i;}
  const Iflt &getr() {return r;}

 private:
  Iflt r,i;
};


#endif

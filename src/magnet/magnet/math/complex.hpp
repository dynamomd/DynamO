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

#include <cmath>

namespace magnet {
  namespace math {
    template <class T>
    class Complex
    {
    public:
      inline Complex():r(0), i(0) {}
      inline Complex(T r2, T i2):r(r2), i(i2) {}
  
      inline Complex &operator+=(const Complex &cc)
      {
	r += cc.r;
	i += cc.i;
	return *this;
      }

      inline Complex &operator-=(const Complex &cc)
      {
	r -= cc.r;
	i -= cc.i;
	return *this;
      }

      inline Complex operator*(const Complex &cc)
      {
	return Complex(r*cc.r-i*cc.i, r*cc.i + i*cc.r);
      }

      inline Complex operator/(const Complex &cc)
      {
	T denom = cc.r * cc.r + cc.i * cc.i;
	return Complex((r * cc.r + i * cc.i)/denom, 
			(i * cc.r - r * cc.i)/denom);
      }

      inline Complex operator*(const T &a)
      {
	return Complex(r*a, i*a);
      }

      inline Complex operator/(const T &a)
      {
	return Complex(r/a, i/a);
      }

      inline bool operator==(const Complex &cc) const
      {
	return ((cc.r == r) && (cc.i == i));
      }

      inline T modulus() const
      {
	T tmp;
	if ((tmp = r*r + i*i) != 0.0)
	  return std::sqrt(tmp);
	else
	  return 0.0;
      }

      inline Complex conjugate()
      {
	return Complex(r,-i);
      }

      inline Complex exponent()
      {
	return Complex(std::exp(r)*std::cos(i), std::exp(r)*std::sin(i));
      }

      inline const T &geti() {return i;}
      inline const T &getr() {return r;}

    private:
      T _r, _i;
    };
  }
}

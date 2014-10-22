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
#include <algorithm>
#include <array>
#include <stddef.h>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <limits>
#include <iomanip>

const size_t NDIM(3);

namespace magnet {
  namespace math {
    template<class T = double, size_t N = 3>
    class NVector : public std::array<T, N> {
      typedef std::array<T, N> Base;
    public:
      NVector() { Base::fill(T());}
      NVector(std::initializer_list<T> _list) {
	if (_list.size() > N)
	  throw std::length_error("initializer list too long");
      
	size_t i = 0;
	auto it = _list.begin();
	for (; it != _list.end(); ++i, ++it)
	  Base::operator[](i) = *it;

	for (; i < N; ++i)
	  Base::operator[](i) = 0.0;
      }

      inline T nrm2() const 
      { 
	T sum(0);
	for (const auto val : *this)
	  sum += val * val;
	return sum;
      }

      inline T nrm() const
      {
	double max = std::abs(Base::operator[](0));
	for (size_t i(1); i < N; ++i)
	  max = std::max(max, std::abs(Base::operator[](i)));
	if (max != 0)
	  {
	    const T invsqrbiggest = 1.0 / (max * max);
	    T sum(0);
	    for (const auto val: *this)
	      sum += invsqrbiggest * val * val;
	    return max * std::sqrt(sum);
	  }
	else
	  return 0;
      }

      inline NVector<T,N> normal() const {
	const T norm = nrm();
	const double inv_nrm = 1.0 / (norm + (norm==0));
	NVector<T,N> retval(*this);
	for (auto& val: retval)
	  val *= inv_nrm;
	return retval;
      }

      inline void normalise() {
	*this = normal();
      }

      inline bool operator==(const NVector<T,N>& ovec) const
      {
	for (size_t i(0); i < N; ++i)
	  if (Base::operator[](i) != ovec[i])
	    return false;
	return true;
      }

      inline bool operator!=(const NVector<T,N>& ovec) const
      {	return !operator==(ovec); }

      template<class P>
      inline NVector<T,N>& operator+=(const P& e)
      { return (*this = *this + e); }
      template<class P>
      inline NVector<T,N>& operator-=(const P& e)
      { return (*this = *this - e); }
      template<class P>
      inline NVector<T,N>& operator*=(const P& e)
      { return (*this = *this * e); }
      template<class P>
      inline NVector<T,N>& operator/=(const P& e)
      { return (*this = *this / e); }

      std::string toString() const
      {
	std::ostringstream os;
	os << std::setprecision(std::numeric_limits<double>::digits10 + 2) << "<";
	for (size_t i(0); i < N-1; ++i)
	  os << Base::operator[](i) << ",";
	os << Base::operator[](N-1) << ">";
	return os.str();
      }
    };
    
    template<class T, size_t N>
    NVector<T,N> operator+(const NVector<T,N>& vec1, const NVector<T,N>& vec2) {
      NVector<T,N> retval(vec1);
      for (size_t i(0); i < N; ++i)
	retval[i] += vec2[i];
      return retval;
    }

    template<class T, size_t N>
    NVector<T,N> operator-(const NVector<T,N>& vec1, const NVector<T,N>& vec2) {
      NVector<T,N> retval(vec1);
      for (size_t i(0); i < N; ++i)
	retval[i] -= vec2[i];
      return retval;
    }

    /*! \brief Outer product (cross) */
    template<class T>
    NVector<T,3> operator^(const NVector<T,3>& vec1, const NVector<T,3>& vec2) {
      return NVector<T,3>{vec1[1] * vec2[2] - vec1[2] * vec2[1],
	  vec1[2] * vec2[0] - vec1[0] * vec2[2],
	  vec1[0] * vec2[1] - vec1[1] * vec2[0]
	  };
    }

    /*! \brief Scalar product (dot) */
    template<class T, size_t N, class R>
    NVector<T,N> operator*(const NVector<T,N>& vec1, const R& val) {
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = vec1[i] * val;
      return retval;
    }
    template<class T, size_t N, class R>
    NVector<T,N> operator*(const R& val, const NVector<T,N>& vec1) {
      return vec1 * val;
    }

    template<class T, size_t N>
    T operator*(const NVector<T,N>& vec1, const NVector<T,N>& vec2) {
      T sum(0);
      for (size_t i(0); i < N; ++i)
	sum += vec1[i] * vec2[i];
      return sum;
    }
    template<class T, size_t N>
    T operator|(const NVector<T,N>& vec1, const NVector<T,N>& vec2) {
      return vec1 * vec2;
    }

    template<class T, size_t N, class R>
    NVector<T,N> operator/(const NVector<T,N>& vec1, const R& val) {
      return vec1 * (R(1)/val);
    }
    template<class T, size_t N, class R>
    NVector<T,N> operator/(const R& val, const NVector<T,N>& vec1) {
      return vec1 * (R(1)/val);
    }

    /*! \brief Unary negative */
    template<class T, size_t N>
    NVector<T,N> operator-(const NVector<T,N>& vec1) {
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = -vec1[i];
      return retval;
    }
    
    template<class T, size_t N>
    inline magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const NVector<T,N>& vec)
    {
      char name[2] = "x";
      for (size_t i(0); i < N; i++)
	{
	  name[0]= 'x'+i; //Write the dimension
	  XML << magnet::xml::attr(name) << vec[i];
	}
      return XML;
    }

    template<class T, size_t N>
    inline NVector<T,N>& operator<<(NVector<T,N>& data, const magnet::xml::Node& XML)
    {
      for (size_t i(0); i < N; i++) 
	{
	  char name[2] = "x";
	  name[0] = 'x'+i; //Write the name
	  data[i] = XML.getAttribute(name).as<T>();
	}

      return data;
    }

    template<class T>
    inline T elementwiseMultiply(const T& A, const T& B)
    { return A * B; }
    
    template<class T, size_t N>
    inline NVector<T,N> elementwiseMultiply(const NVector<T,N>& A, const NVector<T,N>& B)
    { 
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = A[i] * B[i];
      return retval;
    }
    
    template<class T>
    inline T elementwiseMin(const T& A, const T& B)
    { return std::min(A, B); }

    template<class T, size_t N>
    inline NVector<T,N> elementwiseMin(const NVector<T,N>& A, const NVector<T,N>& B)
    { 
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = std::min(A[i], B[i]);
      return retval;
    }
    
    template<class T>
    inline T elementwiseMax(const T& A, const T& B)
    { return std::max(A, B); }

    template<class T, size_t N>
    inline NVector<T,N> elementwiseMax(const NVector<T,N>& A, const NVector<T,N>& B)
    { 
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = std::max(A[i], B[i]);
      return retval;
    }

    typedef NVector<double,3> Vector;
  }
}

namespace coil { typedef ::magnet::math::NVector<double,3> Vector; }
namespace dynamo { typedef ::magnet::math::NVector<double,3> Vector; }

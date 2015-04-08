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
    /*! \brief N-dimensional vector type.
      
      \tparam N The dimensionality of the NVector.

      \tparam T The type of the components of the
      NVector.
     */
    template<class T = double, size_t N = 3>
    class NVector : public std::array<T, N> {
      typedef std::array<T, N> Base;
    public:
      /*! \brief Default constructor.

	This creates a zero vector.
       */
      NVector(T val = T()) { Base::fill(val); }

      /*! \brief Initializer list constructor.

	Allows the simple creation of a vector. For example:
	\code{.cpp}
	NVector<double, 3> vec{0,1,2};
	\endcode	
       */
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

      /*! \brief Conversion copy constructor.
       */
      template<class OT>
      NVector(NVector<OT,N> ov) {
	for (size_t i(0); i < N; ++i)
	  Base::operator[](i) = T(ov[i]);
      }

      /*! \brief Returns the square norm of the NVector.*/
      inline T nrm2() const 
      { 
	T sum(0);
	for (const auto val : *this)
	  sum += val * val;
	return sum;
      }

      /*! \brief Returns the norm of the NVector.*/
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

      /*! \brief Returns the normalised vector in the direction of this NVector.

	NVectors of length zero will return a zero (NVector{0,0,..})
	vector from this function.
       */
      inline NVector<T,N> normal() const {
	const T norm = nrm();
	const double inv_nrm = 1.0 / (norm + (norm==0));
	NVector<T,N> retval(*this);
	for (auto& val: retval)
	  val *= inv_nrm;
	return retval;
      }

      /*! \brief Normalise this NVector.

	For a zero NVector, this will not alter the vector.
       */
      inline void normalise() {
	*this = normal();
      }

      /*! \brief Comparison operator.*/
      inline bool operator==(const NVector<T,N>& ovec) const
      {
	for (size_t i(0); i < N; ++i)
	  if (Base::operator[](i) != ovec[i])
	    return false;
	return true;
      }

      /*! \brief Comparison operator.*/
      inline bool operator!=(const NVector<T,N>& ovec) const
      {	return !operator==(ovec); }

      /*! \name Modify-assign operators

	These are for convenience only and are not optimised.
	\{
       */
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
      /*! \} */

      /*! \brief Create a human-readable representation of this NVector.*/
      std::string toString() const
      {
	std::ostringstream os;
	os << std::setprecision(std::numeric_limits<double>::digits10 + 2) << "Vector{";
	for (size_t i(0); i < N-1; ++i)
	  os << Base::operator[](i) << ",";
	os << Base::operator[](N-1) << "}";
	return os.str();
      }
    };

    /*! \relates NVector
      \name Vector Arithmetic 
      \{
     */
    /*! \brief Addition of two NVector types.  */
    template<class T1, class T2, size_t N>
    NVector<decltype(T1()+T2()), N> operator+(const NVector<T1,N>& vec1, const NVector<T2,N>& vec2) {
      NVector<decltype(T1()+T2()), N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = vec1[i] + vec2[i];
      return retval;
    }

    /*! \brief Subtraction of two NVector types.  */
    template<class T1, class T2, size_t N>
    NVector<decltype(T1()-T2()),N> operator-(const NVector<T1,N>& vec1, const NVector<T2,N>& vec2) {
      NVector<decltype(T1()-T2()),N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = vec1[i] - vec2[i];
      return retval;
    }

    /*! \brief Outer (cross) product  */
    template<class T>
    NVector<T,3> operator^(const NVector<T,3>& vec1, const NVector<T,3>& vec2) {
      return NVector<T,3>{vec1[1] * vec2[2] - vec1[2] * vec2[1],
	  vec1[2] * vec2[0] - vec1[0] * vec2[2],
	  vec1[0] * vec2[1] - vec1[1] * vec2[0]
	  };
    }

    /*! \brief Multiplication of a scalar and an NVector.  */
    template<class T, size_t N, class R,
	     typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>
    NVector<decltype(T()*R()),N> operator*(const NVector<T,N>& vec1, const R& val) {
      NVector<decltype(T()*R()),N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = vec1[i] * val;
      return retval;
    }

    /*! \brief Multiplication of an NVector and a scalar.  */
    template<class T, size_t N, class R,
	     typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>
    NVector<decltype(R()*T()),N> operator*(const R& val, const NVector<T,N>& vec1) {
      NVector<decltype(R()*T()),N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = val * vec1[i];
      return retval;
    }

    /*! \brief Scalar (dot) product */
    template<class T1, class T2, size_t N>
    decltype(T1()*T2()) operator*(const NVector<T1,N>& vec1, const NVector<T2,N>& vec2) {
      decltype(T1()*T2()) sum(0);
      for (size_t i(0); i < N; ++i)
	sum += vec1[i] * vec2[i];
      return sum;
    }

    /*! \brief Scalar (dot) product */
    template<class T1, class T2, size_t N>
    decltype(NVector<T1,N>()* NVector<T2,N>()) operator|(const NVector<T1,N>& vec1, const NVector<T2,N>& vec2) {
      return vec1 * vec2;
    }

    /*! \brief Division of an NVector by a scalar.  */
    template<class T, size_t N, class R,
	     typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>
    NVector<decltype(T()/R()),N> operator/(const NVector<T,N>& vec1, const R& val) {
      NVector<decltype(T()/R()),N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = vec1[i] / val;
      return retval;
    }

    /*! \brief Unary negative operator. */
    template<class T, size_t N>
    NVector<T,N> operator-(const NVector<T,N>& vec1) {
      NVector<T,N> retval;
      for (size_t i(0); i < N; ++i)
	retval[i] = -vec1[i];
      return retval;
    }
    
    /*! \} */

    /*! \relates NVector
      \name Vector input/output operators
      \{
     */

    /*! \brief XML output operator. */
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

    /*! \brief XML input operator. */
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

    /*! \brief Output a human readable representation to an output stream. */
    template<class Real, size_t N>
    inline std::ostream& operator<<(std::ostream& os, const NVector<Real, N>& vec) {
      os << vec.toString();
      return os;
    }

    /*! \} */

    /*! \relates NVector
      \name Elementwise operations on NVectors
      \{
     */

    template<class T1, class T2>
    inline decltype(T1()*T2()) elementwiseMultiply(const T1& A, const T2& B)
    { return A * B; }
    
    template<class T1, class T2, size_t N>
    inline NVector<decltype(T1()*T2()),N> elementwiseMultiply(const NVector<T1,N>& A, const NVector<T2,N>& B)
    { 
      NVector<decltype(T1()*T2()),N> retval;
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

    /*! \} */    

    typedef NVector<double,3> Vector;
  }
}

namespace coil { typedef ::magnet::math::NVector<double,3> Vector; }
namespace dynamo { typedef ::magnet::math::NVector<double,3> Vector; }

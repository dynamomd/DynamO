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

#ifndef VECTOR_H
#define VECTOR_H

#include <cmath> //for length
#include "../base/is_exception.hpp"
#include "../base/constants.hpp" //For the dimensions

class XMLNode;

namespace xmlw 
{
  class XmlStream;
}

template <typename T = Iflt> 
class CVector
{
 public:
  T data[NDIM];

  inline CVector() {}
    
  inline CVector (const T& fillItem)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] = fillItem;
  }
  
  inline CVector<T> operator+ (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] + v2.data[i];
    
    return tmpVec;
  }
  
  inline CVector<T>& operator+= (const CVector < T > &v2)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] += v2.data[i];
    
    return *this;
  }
  
  inline CVector < T > operator- (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] - v2.data[i];
    
    return tmpVec;
  }
  
  inline CVector<T>& operator-= (const CVector < T > &v2)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] -= v2.data[i];
    
    return *this;
  }
  
  inline CVector<T> operator* (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] * v2.data[i];
    
    return tmpVec;
  }
  
  inline CVector<T>& operator*= (const CVector < T > &v2)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] *= v2.data[i];
    
    return *this;
  }
  
  inline CVector<T> operator* (const T &val) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] * val;
    
    return tmpVec;
  }
  
  friend inline CVector<T> operator* (const T &val,
				      const CVector<T>& vec)
  { return vec * val; }
  
  inline CVector<T>& operator*= (const T &val)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] *= val;
    
    return *this;
  }
  
  inline CVector<T> operator/ (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] / v2.data[i];
    
    return tmpVec;
  }
  
  inline CVector < T > operator/ (const T &val) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] / val;
    
    return tmpVec;
  }
  
  inline CVector<T> operator/= (const T &val)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] /= val;
    
    return *this;
  }
  
  inline CVector<T> Cross(const CVector<T> &val)
  {
#ifdef DYNAMO_DEBUG
    if (NDIM != 3) 
      I_throw() << "Cross product defined only in 3D";
#endif
    
    CVector < T > tmpVec;
    tmpVec[0] = data[1] * val.data[2] -  data[2] * val.data[1];
    tmpVec[1] = data[2] * val.data[0] -  data[0] * val.data[2];
    tmpVec[2] = data[0] * val.data[1] -  data[1] * val.data[0];
    
    return tmpVec;
  }
  
  // The dot product operator
  inline T operator% (const CVector < T > &v2) const
  {
    T tmpVal = data[0] * v2.data[0];
    for (int iDim = 1; iDim < NDIM; ++iDim)
      tmpVal += data[iDim] * v2.data[iDim];
    
    return tmpVal;
  }
  
  inline CVector < T > unitVector() const
  {
    return (*this) / length();
  }
  
  inline T square() const
  {
    T tmpVal(data[0] * data[0]);
    
    for (int iDim = 1; iDim < NDIM; ++iDim)
      tmpVal += data[iDim] * data[iDim];
    
    return tmpVal;
  }
  
  inline T length() const
  { return std::sqrt(square()); }
  
  inline CVector < T > operator-()const
  { 
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = -data[i];
    
    return tmpVec;
  }
  
  inline T & operator[] (const int &indx)
  {
    return data[indx];
  }
  
  inline const T & operator[] (const int &indx) const
  {
    return data[indx];
  }
  
  inline CVector<CVector<T> > dyad(const CVector<T>& rightVec) const
  {    
    CVector<CVector<T> > tmpvec(0);

    for (int iDim = 0; iDim < NDIM; iDim++)
      for (int jDim = 0; jDim < NDIM; jDim++)
	tmpvec[iDim][jDim] = data[iDim] * rightVec.data[jDim];	\
    
    return tmpvec;
  }
  
  void operator<<(const XMLNode&);
  
  template<class S>
  friend xmlw::XmlStream& 
  operator<<(xmlw::XmlStream&, const CVector<S>&);
  
  template<class S>
  friend xmlw::XmlStream& 
  operator<<(xmlw::XmlStream&, const CVector<CVector<S> >&);

  template<class A>
  friend CVector<A>
  convertVec(const CVector<T>& old)
  {
    CVector<A> tmp;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      tmp[iDim] = static_cast<A>(old[iDim]);

    return tmp;
  }

};

#endif

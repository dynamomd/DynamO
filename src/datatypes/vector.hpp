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
/*! \file vector.hpp
 * Contains the non XML functions and definitions of the CVector class.
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

/*! \brief The simple NDIM vector class.
 *
 * This class is used for vector storage and operations. It seems
 * to optimise to zero/low overhead as it has resisted all my attempts to
 * optimise it further. I've currently tried templated/preprocessor unrolling
 * and I failed at using liboil to take advantage of the systems MMX
 * extensions etc.
 */
template <typename T = Iflt> 
class CVector
{
 public:
  /*! \brief The data of the vector.
   */
  T data[NDIM];

  /*! \brief The default initialisor for speed. */
  inline CVector() {}
    
  /*! \brief Initialise all elements of the vector to a certain value. */
  inline explicit CVector(const T& fillItem)
  {
    for (size_t i(0); i < NDIM; ++i)
      data[i] = fillItem;
  }
  
  /*! \brief Adds each dimension of the vectors seperately. */
  inline CVector<T> operator+ (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] + v2.data[i];
    
    return tmpVec;
  }
  
  /*! \brief Adds each dimension of the vectors seperately. */
  inline CVector<T>& operator+= (const CVector < T > &v2)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] += v2.data[i];
    
    return *this;
  }
  
  /*! \brief Subtracts each dimension of the vectors seperately. */
  inline CVector < T > operator- (const CVector < T > &v2) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] - v2.data[i];
    
    return tmpVec;
  }
  
  /*! \brief Subtracts each dimension of the vectors seperately. */
  inline CVector<T>& operator-= (const CVector < T > &v2)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] -= v2.data[i];
    
    return *this;
  }
  
  /*! \brief Scales the vector by a factor */
  inline CVector<T> operator* (const T &val) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] * val;
    
    return tmpVec;
  }
  
  /*! \brief Scales the vector by a factor */
  friend inline CVector<T> operator* (const T &val,
				      const CVector<T>& vec)
  { return vec * val; }
  
  /*! \brief Scales the vector by a factor */
  inline CVector<T>& operator*= (const T &val)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] *= val;
    
    return *this;
  }
  
  /*! \brief Divides the vector by a factor */
  inline CVector < T > operator/ (const T &val) const
  {
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = data[i] / val;
    
    return tmpVec;
  }
  
  /*! \brief Divides the vector by a factor */
  inline CVector<T> operator/= (const T &val)
  {
    for (size_t i = 0; i < NDIM; ++i)
      data[i] /= val;
    
    return *this;
  }
  
  /*! \brief Only defined for 3d but outputs the cross product. */  
  inline CVector<T> Cross(const CVector<T> &val)
  {
#ifdef DYNAMO_DEBUG
    if (NDIM != 3) 
      D_throw() << "Cross product defined only in 3D";
#endif
    
    CVector < T > tmpVec;
    tmpVec[0] = data[1] * val.data[2] -  data[2] * val.data[1];
    tmpVec[1] = data[2] * val.data[0] -  data[0] * val.data[2];
    tmpVec[2] = data[0] * val.data[1] -  data[1] * val.data[0];
    
    return tmpVec;
  }
  
  /*! \brief The dot product operator. */
  inline T operator% (const CVector < T > &v2) const
  {
    T tmpVal = data[0] * v2.data[0];
    for (int iDim = 1; iDim < NDIM; ++iDim)
      tmpVal += data[iDim] * v2.data[iDim];
    
    return tmpVal;
  }
  
  /*! \brief Returns a unit CVector of the current CVector*/
  inline CVector < T > unitVector() const
  {
    return (*this) / length();
  }
  
  /*! \brief Returns the value of the vector dotted with itself.*/
  inline T square() const
  {
    T tmpVal(data[0] * data[0]);
    
    for (int iDim = 1; iDim < NDIM; ++iDim)
      tmpVal += data[iDim] * data[iDim];
    
    return tmpVal;
  }
  
  /*! \brief Returns the scalar length of the vector.*/
  inline T length() const
  { return std::sqrt(square()); }
  
  /*! \brief Flips the CVector's direction. */
  inline CVector < T > operator-()const
  { 
    CVector < T > tmpVec;
    
    for (size_t i = 0; i < NDIM; ++i)
      tmpVec.data[i] = -data[i];
    
    return tmpVec;
  }
  
  /*! \brief Non-const accessor to a dimension of the vector */
  inline T & operator[] (const size_t &indx)
  {
#ifdef DYNAMO_DEBUG
    if (indx >= NDIM) D_throw() << "CVector out of bounds error";
#endif

    return data[indx];
  }
  
  /*! \brief const accessor to a dimension of the vector */
  inline const T & operator[] (const size_t &indx) const
  {
#ifdef DYNAMO_DEBUG
    if (indx >= NDIM) D_throw() << "CVector out of bounds error";
#endif

    return data[indx];
  }
  
  /*! \brief Returns the dyadic product of two vectors.  */  
  inline CVector<CVector<T> > dyad(const CVector<T>& rightVec) const
  {    
    CVector<CVector<T> > tmpvec(CVector<>(0));

    for (int iDim = 0; iDim < NDIM; iDim++)
      for (int jDim = 0; jDim < NDIM; jDim++)
	tmpvec[iDim][jDim] = data[iDim] * rightVec.data[jDim];	\
    
    return tmpvec;
  }
  
  /*! \brief Will load a CVector from an XMLNode.*/
  void operator<<(const XMLNode&);
  
  /*! \brief A helper function to write out a CVector to an XMLStream*/
  template<class S>
  friend xmlw::XmlStream& 
  operator<<(xmlw::XmlStream&, const CVector<S>&);
  
  /*! \brief A helper function to write out a CVector matrix to an XMLStream*/
  template<class S>
  friend xmlw::XmlStream& 
  operator<<(xmlw::XmlStream&, const CVector<CVector<S> >&);

  /*! \brief A helper function to change the stored type of a vector. */
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

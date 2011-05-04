/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
/*! \file vector.xml.hpp
 *
 * Contains the definition of the XML functions for a CVector
 */

#pragma once

#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/math/vector.hpp>
#include "vector.hpp"

template<class T>
void 
CVector<T>::operator<<(const magnet::xml::Node& XML)
{
  for (int iDim = 0; iDim < NDIM; iDim++) 
    {
      char name[2] = "x";
      name[0] = 'x' + iDim; //Write the name
      if (!XML.getAttribute(name).valid())
	name[0] = '0'+iDim;
      
      try {
	data[iDim] = XML.getAttribute(name).as<T>();
      }
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in CVector";
	}
    }
}

template<class T>
inline
xml::XmlStream& 
operator<<(xml::XmlStream& XML, const CVector<T> &vec)
{
  char name[2] = "x";
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xml::attr(name) << vec[iDim];
    }
  
  return XML;
}

template<class T>
inline
xml::XmlStream& 
operator<<(xml::XmlStream& XML, const CVector<CVector<T> > &vec)
{
  char name[2] = "x";
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xml::tag(name) << vec[iDim]
	  << xml::endtag(name);
    }
  
  return XML;
}

// vectors
template<class A, int B, class C>
inline xml::XmlStream& operator<<(xml::XmlStream& XML, 
				    const VectorExpression<A,B,C> & t )
{
  char name[2] = "x";
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xml::attr(name) << t(iDim);
    }
  
  return XML;
}

inline
VectorExpression<>& 
operator<<(VectorExpression<>& data, const magnet::xml::Node& XML)
{
  for (size_t iDim = 0; iDim < NDIM; iDim++) 
    {
      char name[2] = "x";
      name[0] = 'x' + iDim; //Write the name
      if (!XML.getAttribute(name))
	name[0] = '0'+iDim;
      
      try {
	data[iDim] = XML.getAttribute(name).as<double>();
      }
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in CVector";
	}
    }

  return data;
}

// vectors
template<class A, int B, class C>
inline xml::XmlStream& operator<<(xml::XmlStream& XML, 
				   const MatrixExpression<A,B,C> & t )
{
  char name[2] = "x";

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      name[0] = 'x'+iDim;
      XML << xml::tag(name);
      
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  char name2[2] = "x";
	  name2[0] = 'x'+jDim;
	  XML << xml::attr(name2) << t(iDim,jDim);
	}

      XML << xml::endtag(name);
    }
  
  return XML;
}

#ifdef MATRIX_HEADER
inline
MatrixExpression<>& 
operator<<(MatrixExpression<>& data, const magnet::xml::Node& XML)
{
  char name[2] = "x";
  
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      name[0] = 'x'+iDim;

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  char name2[2] = "x";
	  name2[0] = 'x'+jDim;

	  try {
	    data(iDim,jDim) = XML.getNode(name).getAttribute(name2).as<double>();
	  }
	  catch (boost::bad_lexical_cast &)
	    {
	      M_throw() << "Failed a lexical cast";
	    }
	}
    }

  return data;
}
#endif

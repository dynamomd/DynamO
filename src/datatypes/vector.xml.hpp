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
/*! \file vector.xml.hpp
 *
 * Contains the definition of the XML functions for a CVector
 */

#ifndef VECTOR_XML_H
#define VECTOR_XML_H

#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include "../base/constants.hpp"
#include <boost/lexical_cast.hpp> //For xml Parsing

template<class T>
void 
CVector<T>::operator<<(const XMLNode &XML)
{
  for (int iDim = 0; iDim < NDIM; iDim++) 
    {
      char name[2] = "x";
      name[0] = 'x' + iDim; //Write the name
      if (!XML.isAttributeSet(name))
	name[0] = '0'+iDim;
      
      try {
	data[iDim] = boost::lexical_cast<T>(XML.getAttribute(name));
      }
      catch (boost::bad_lexical_cast &)
	{
	  D_throw() << "Failed a lexical cast in CVector";
	}
    }
}

template<class T>
inline
xmlw::XmlStream& 
operator<<(xmlw::XmlStream& XML, const CVector<T> &vec)
{
  char name[2] = "x";
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xmlw::attr(name) << vec[iDim];
    }
  
  return XML;
}

template<class T>
inline
xmlw::XmlStream& 
operator<<(xmlw::XmlStream& XML, const CVector<CVector<T> > &vec)
{
  char name[2] = "x";
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xmlw::tag(name) << vec[iDim]
	  << xmlw::endtag(name);
    }
  
  return XML;
}

// vectors
template<class A, int B, class C>
inline xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
				    const VectorExpression<A,B,C> & t )
{
  char name[2] = "x";
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      name[0]= 'x'+iDim; //Write the dimension
      XML << xmlw::attr(name) << t(iDim);
    }
  
  return XML;
}

inline
VectorExpression<>& 
operator<<(VectorExpression<>& data, const XMLNode &XML)
{
  for (size_t iDim = 0; iDim < NDIM; iDim++) 
    {
      char name[2] = "x";
      name[0] = 'x' + iDim; //Write the name
      if (!XML.isAttributeSet(name))
	name[0] = '0'+iDim;
      
      try {
	data[iDim] = boost::lexical_cast<DOUBLE>(XML.getAttribute(name));
      }
      catch (boost::bad_lexical_cast &)
	{
	  D_throw() << "Failed a lexical cast in CVector";
	}
    }

  return data;
}

// vectors
template<class A, int B, class C>
inline xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
				   const MatrixExpression<A,B,C> & t )
{
  char name[2] = "x";

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      name[0] = 'x'+iDim;
      XML << xmlw::tag(name);
      
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  char name2[2] = "x";
	  name2[0] = 'x'+jDim;
	  XML << xmlw::attr(name2) << t(iDim,jDim);
	}

      XML << xmlw::endtag(name);
    }
  
  return XML;
}

inline
MatrixExpression<>& 
operator<<(MatrixExpression<>& data, const XMLNode &XML)
{
  char name[2] = "x";
  
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      name[0] = 'x'+iDim;

      XMLNode BrowseNode = XML.getChildNode(name);

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  char name2[2] = "x";
	  name2[0] = 'x'+jDim;

	  try {
	    data(iDim,jDim) = boost::lexical_cast<DOUBLE>(BrowseNode.getAttribute(name2));
	  }
	  catch (boost::bad_lexical_cast &)
	    {
	      D_throw() << "Failed a lexical cast in CVector";
	    }
	}
    }

  return data;
}

#endif

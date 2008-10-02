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

#ifndef BCShapes_H
#define BCShapes_H

#include "BC.hpp"
#include "../../base/is_simdata.hpp"
#include <cmath> 
#include <boost/preprocessor.hpp>

class CRectBC: public CBC 
{
protected:
  //Not virtual as here ends the rounding definition
  inline void rounding(CVector<>& pos) const
  {

#ifdef DYNAMO_double_precsision
# define BOOST_PP_LOCAL_LIMITS (0, NDIM - 1)
# define BOOST_PP_LOCAL_MACRO(n)			\
    pos.data[n] -= Sim->aspectRatio[n] *	\
      rint (pos.data[n]/Sim->aspectRatio[n]);
    /**/
    
# include BOOST_PP_LOCAL_ITERATE()
#else
# define BOOST_PP_LOCAL_LIMITS (0, NDIM - 1)
# define BOOST_PP_LOCAL_MACRO(n)			\
    pos.data[n] -= Sim->aspectRatio[n] *	\
      rintf (pos.data[n]/Sim->aspectRatio[n]);
    /**/
    
# include BOOST_PP_LOCAL_ITERATE()
#endif
  }
};

class CSqBC: public CBC
{  
protected:
  //Not virtual as here ends the rounding definition
  inline void rounding(CVector<>& pos) const
  {
#ifdef DYNAMO_double_precsision
# define BOOST_PP_LOCAL_LIMITS (0, NDIM - 1)
# define BOOST_PP_LOCAL_MACRO(n)			\
      pos.data[n] -= rint (pos.data[n]);

# include BOOST_PP_LOCAL_ITERATE()
#else
# define BOOST_PP_LOCAL_LIMITS (0, NDIM - 1)
# define BOOST_PP_LOCAL_MACRO(n)			\
      pos.data[n] -= rintf (pos.data[n]);

# include BOOST_PP_LOCAL_ITERATE()
#endif
  }
};


#endif

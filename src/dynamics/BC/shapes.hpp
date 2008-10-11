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
#include "../../base/is_exception.hpp"

#ifdef DYNAMO_double_precsision
# define rintfunc rint
#else
# define rintfunc rintf
#endif

/*! \brief Just a simple class that implements rectangular system rounding.
 *
 * The definition between square and rectangular is only for
 * optimisation purposes.
 */
class CRectBC: public CBC 
{
public:
  CRectBC(const DYNAMO::SimData* const &SD, const char *aName, 
      const char *aColor):
    CBC(SD, aName, aColor)
  {
  }

protected:
  //Not virtual as here ends the rounding definition
  inline void rounding(CVector<>& pos) const
  {
    for (size_t n = 0; n < NDIM; ++n)
      pos.data[n] -= Sim->aspectRatio[n] *
	rintfunc (pos.data[n]/Sim->aspectRatio[n]);    
  }
};

/*! \brief Just a simple class that implements square system rounding.
 *
 * The definition between square and rectangular is only for
 * optimisation purposes.
 */
class CSqBC: public CBC
{  
public:
  CSqBC(const DYNAMO::SimData* const &SD, const char *aName, 
	const char *aColor):
    CBC(SD, aName, aColor)
  {
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if (Sim->aspectRatio[iDim] != 0)
	D_throw() << "The simulation aspect ratio is not unity for the use of "
	  "square PBC's";
  }

protected:
  //Not virtual as here ends the rounding definition
  inline void rounding(CVector<>& pos) const
  {
    for (size_t n = 0; n < NDIM; ++n)
      pos.data[n] -= rintfunc (pos.data[n]);
  }
};


#endif

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

#pragma once
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <cmath>

struct CUringRod: public CUCell
{
  CUringRod(size_t pcl, Iflt WL, CUCell* nextCell):
    CUCell(nextCell),
    pairchainlength(pcl),
    walklength(WL)
  {
    if (pcl == 0) I_throw() << "Cant have zero chain length";
  }

  size_t pairchainlength;  
  Iflt walklength;
  
  virtual std::vector<CVector<> > placeObjects(const CVector<>& centre)
  {
    std::vector<CVector<> > localsites;
        
    for (size_t iStep = 0; iStep < pairchainlength; ++iStep)
      { 
	CVector<> tmp(0);
	tmp[0] = -0.5 * walklength;
	tmp[1] = walklength * ( iStep - 0.5 * (pairchainlength-1));

	localsites.push_back(tmp + centre);
      }

    for (int iStep = pairchainlength; iStep != 0;)
      { 
	CVector<> tmp(0);

	--iStep;
	tmp[0] = 0.5 * walklength;
	tmp[1] = walklength * ( iStep - 0.5 * (pairchainlength-1));

	localsites.push_back(tmp + centre);
      }
  
    std::vector<CVector<> > retval;
    BOOST_FOREACH(const CVector<>& vec, localsites)
      BOOST_FOREACH(const CVector<>& vec2, uc->placeObjects(vec))
        retval.push_back(vec2);

    return retval;    
  }
};

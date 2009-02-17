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

struct CUHelix: public CUCell
{
  CUHelix(long CL, long RL, Iflt WL, Iflt D, CUCell* nextCell):
    CUCell(nextCell),
    chainlength(CL),
    ringlength(RL),
    walklength(WL),
    diameter(D)
  {}

  long chainlength;
  long ringlength;
  Iflt walklength;
  Iflt diameter;
  
  virtual std::vector<CVector<> > placeObjects(const CVector<>& centre)
  {
    Iflt a = diameter * (0.5 / PI);
    Iflt sigstep = 2.0 * PI / ringlength;
    Iflt zcentre = a * (chainlength - 1) * sigstep * a;
    Iflt radius = 0.5 * std::sqrt(walklength * walklength - std::pow(a / ringlength,2)) / std::sin(PI / ringlength);

    std::vector<CVector<> > localsites;
    
    CVector<> tmp;
    
    for (int iStep = 0; iStep < chainlength; ++iStep)
      { 
	tmp[0] = radius * std::cos(sigstep * iStep);
	tmp[1] = radius * std::sin(sigstep * iStep);
	tmp[2] = a * sigstep * iStep - zcentre;

	localsites.push_back(tmp + centre);
      }
  
    std::vector<CVector<> > retval;
    BOOST_FOREACH(const CVector<>& vec, localsites)
      BOOST_FOREACH(const CVector<>& vec2, uc->placeObjects(vec))
        retval.push_back(vec2);

    return retval;    
  }
};

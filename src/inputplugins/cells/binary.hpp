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

struct CUBinary: public CUCell
{
  CUBinary(size_t x, CUCell* nextCell1, CUCell* nextCell2):
    CUCell(nextCell1),
    uc2(nextCell2),
    count(0),
    countA(x)
  {}

  boost::scoped_ptr<CUCell> uc2;
  size_t count;
  const size_t countA;
  
  virtual std::vector<Vector> placeObjects(const Vector & centre)
  {
    if (count < countA)
      {
	++count;
	return uc->placeObjects(centre);
      }
    else
	return uc2->placeObjects(centre);
  }

};

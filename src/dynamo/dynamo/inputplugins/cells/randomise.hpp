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

#pragma once
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <algorithm>

struct CURandomise: public CUCell
{
  CURandomise(CUCell* nextCell):
    CUCell(nextCell)
  {}

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    //Must be placed at zero for the mirroring to work correctly
    std::vector<Vector  > retval(uc->placeObjects(Vector (centre)));
    
    std::random_shuffle(retval.begin(), retval.end());
    
    return retval;    
  }
};

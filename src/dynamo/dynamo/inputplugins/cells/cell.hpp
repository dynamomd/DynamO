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
#include <vector>
#include "../../datatypes/vector.hpp"
#include <boost/scoped_ptr.hpp>

class CUCell
{
public:
  CUCell(CUCell* nP): uc(nP) {}

  virtual ~CUCell() {}

  virtual void initialise() { uc->initialise(); }

  virtual std::vector<Vector  > placeObjects(const Vector & ) = 0;  
  
  boost::scoped_ptr<CUCell> uc;
};

//A simple terminator, used to place a particle at this point
class CUParticle: public CUCell
{
public:
  CUParticle(): CUCell(NULL) {}

  //Terminate initialisation
  virtual void initialise() {}

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    std::vector<Vector  > retval;
    retval.push_back(centre);
    return retval;
  }
};

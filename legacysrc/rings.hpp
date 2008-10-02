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

#ifndef CIRings_H
#define CIRings_H

#include "inputplugin.hpp"
#include "../datatypes/vector.hpp"
#include <list>

class CIRings : public CInputPlugin
{
 public:
  CIRings(Iflt, CVector<long>, long, DYNAMO::SimData*); 

  virtual void setSimType(unsigned int);

  void initialise();

 protected:
  
  Iflt density, volume;
  CVector<long> cells;
  CVector<> aspectRatio;
  long maxdim, ncells;
  unsigned int chainlength;

  Iflt latticeWidth;
  Iflt systemWidth;
  Iflt siteAngle;
  Iflt siteRadius;
  Iflt atomDiam;
  Iflt bondRadius;
  Iflt bondwidth;

  std::list<CVector<> > sitelists;
};
#endif

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
#include "../base/is_base.hpp"
#include "fuzzy_array.hpp"
#include "../datatypes/vector.hpp"
#include "../base/is_simdata.hpp"

#define NBins 32

template<class T>
class CFieldArray: public dynamo::SimBase_const
{
 public:
  CFieldArray(const dynamo::SimData* tmp):
    SimBase_const(tmp,"FieldArray",IC_cyan),
    Field(1.0/NBins, -0.5, NBins)
    {};
  
  template<class S>
    CFuzzyArray2<CFuzzyArray2<T> > & operator[](const S x)
    { return Field[x]; }

  T& operator[](const Vector  &cv)
    { 
      return Field
	[cv.data[0]/Sim->primaryCellSize[0]]
	[cv.data[1]/Sim->primaryCellSize[1]]
	[cv.data[2]/Sim->primaryCellSize[2]]; 
    }
    
  long getnBins() const
  { return Field.data.size(); }
  
  T getAverage() 
  {
    T sum = 0;
    for (long z = 0; z < NBins; z++)
      for (long y = 0; y < NBins; y++)
	for (long x = 0; x < NBins; x++)
	  sum += Field[x][y][z];
    sum /= NBins*NBins*NBins;
    return sum;
  }

  CFuzzyArray2<CFuzzyArray2<CFuzzyArray2<T> > > Field;
};  

/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/inputplugins/cells/cell.hpp>
#include <random>

namespace dynamo {
  struct CURandWalk: public UCell
  {
    CURandWalk(long CL, double WL, double D, UCell* nextCell):
      UCell(nextCell),
      chainlength(CL),
      walklength(WL),
      diameter(D),
      _rng(std::random_device()())
    {}

    long chainlength;
    double walklength;
    double diameter;
  
    std::mt19937 _rng;

    Vector getRandVec()
    {
      //See http://mathworld.wolfram.com/SpherePointPicking.html
      std::normal_distribution<> normdist(0.0, (1.0 / sqrt(double(NDIM))));
    
      Vector  tmpVec;
      for (size_t iDim = 0; iDim < NDIM; iDim++)
	tmpVec[iDim] = normdist(_rng);
    
      return tmpVec;
    }

    virtual std::vector<Vector> placeObjects(const Vector & centre)
    {
      std::vector<Vector> localsites;
    
      Vector start{0,0,0}, tmp{0,0,0};
    
      for (int iStep = 0; iStep < chainlength; ++iStep)
	{      
	  bool test = true;
	
	  while(test)
	    {
	      test = false;
	    
	      {
		Vector tmp2(getRandVec());
		tmp2 /= tmp2.nrm();
	      
		tmp = start + tmp2 * walklength;
	      }
	    
	      for (const Vector & vec : localsites)
		if ((vec - tmp).nrm() <= diameter)
		  test = true;
	    }
		
	  localsites.push_back(start);
	
	  start = tmp;
	}
  
      //Centre the chain in the unit cell
      {
	Vector offset{0,0,0};
      
	for (const Vector & vec : localsites)
	  offset += vec;
      
	offset /= localsites.size();
      
	// move to the centre and offset
	for (Vector & vec : localsites)
	  vec -= offset + centre;
      }

      std::vector<Vector> retval;
      for (const Vector& vec : localsites)
	{
	  const std::vector<Vector>& newsites = uc->placeObjects(vec);
	  retval.insert(retval.end(), newsites.begin(), newsites.end());
	}

      return retval;    
    }
  };
}

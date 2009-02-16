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
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

struct CURandWalk: public CUCell
{
  CURandWalk(long CL, Iflt WL, Iflt D, CUCell* nextCell):
    CUCell(nextCell),
    chainlength(CL),
    walklength(WL),
    diameter(D),
    ranGenerator(static_cast<unsigned>(std::time(0)))
  {}

  long chainlength;
  Iflt walklength;
  Iflt diameter;
  
  boost::mt19937 ranGenerator;

  CVector<> getRandVelVec()
  {
    //See http://mathworld.wolfram.com/SpherePointPicking.html
    boost::normal_distribution<Iflt> normdist(0.0, (1.0 / sqrt(NDIM)));
    
    boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution<Iflt> >
      normal_sampler(ranGenerator, normdist);
    
    CVector<> tmpVec;
    for (int iDim = 0; iDim < NDIM; iDim++)
      tmpVec[iDim] = normal_sampler();
    
    return tmpVec;
  }

  virtual std::vector<CVector<> > placeObjects(const CVector<>& centre)
  {
    std::vector<CVector<> > localsites;
    
    CVector<> start(0.0), tmp(0.0);
    
    for (int iStep = 0; iStep < chainlength; ++iStep)
      {      
	bool test = true;
	
	while(test)
	  {
	    test = false;
	    
	    tmp = start + getRandVelVec().unitVector() * walklength;
	    
	    BOOST_FOREACH(const CVector<>& vec, localsites)
	      if ((vec - tmp).length() <= diameter)
		test = true;
	  }
		
	localsites.push_back(start);
	
	start = tmp;
      }
  
    //Centre the chain in the unit cell
    {
      CVector<> offset(0.0);
      
      BOOST_FOREACH(const CVector<>& vec, localsites)
	offset += vec;
      
      offset /= localsites.size();
      
      // move to the centre and offset
      BOOST_FOREACH(CVector<>& vec, localsites)
	vec -= offset + centre;
    }

    std::vector<CVector<> > retval;
    BOOST_FOREACH(const CVector<>& vec, localsites)
      BOOST_FOREACH(const CVector<>& vec2, uc->placeObjects(vec))
        retval.push_back(vec2);

    return retval;    
  }
};

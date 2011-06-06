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
#include <magnet/math/vector.hpp>
#include <magnet/math/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <cmath>

struct CUCylinder: public CUCell
{
  CUCylinder(double partD, double cylD, Vector naxis,
	     boost::uniform_01<dynamo::baseRNG, double>& rng, 
	     CUCell* nextCell):
    CUCell(nextCell),
    diameter(cylD),
    minSpacing(partD),
    axis(naxis),
    uniformRng(rng)
  {}

  double diameter;
  double minSpacing;
  Vector axis;
  boost::uniform_01<dynamo::baseRNG, double>& uniformRng;

  virtual std::vector<Vector> placeObjects(const Vector & centre)
  {
    Vector startpoint = centre - 0.5 * axis;
    
    Vector perpvector;
    for (size_t id(0); id < NDIM; ++id)
      perpvector[id] = uniformRng();

    perpvector -=  axis * ((perpvector | axis) / axis.nrm2());

    perpvector *= 0.5 * diameter / perpvector.nrm();
    
    Vector unitAxis = axis / axis.nrm();

    std::vector<Vector> localsites;
    
    size_t nperrad = M_PI * diameter / minSpacing;
    size_t nrings = axis.nrm() / minSpacing;
    double arcsize = 2.0 * M_PI / nperrad;
    double ringstep = axis.nrm() / nrings;

    for (size_t iStep = 0; iStep < nrings; ++iStep)
      {
	for (size_t jStep = 0; jStep < nperrad; ++jStep)
	  { 
	    Vector pos = startpoint 
	      + (Rodrigues(unitAxis * jStep * arcsize) * perpvector);
	      
	    localsites.push_back(pos);
	  }
	startpoint += unitAxis * ringstep;
      }
  
    std::vector<Vector> retval;

    BOOST_FOREACH(const Vector& vec, localsites)
      BOOST_FOREACH(const Vector& vec2, uc->placeObjects(vec))
        retval.push_back(vec2);
    
    return retval;    
  }
};

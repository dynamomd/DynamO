/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <cmath>
#include <limits>
#include <magnet/math/quadratic.hpp>

namespace magnet {
  namespace math {
    //Please read  http://linus.it.uts.edu.au/~don/pubs/solving.html
    //For solving cubics like x^3 + p * x^2 + q * x + r == 0
    //This function always returns root1 >= root2 >= root3!
    inline size_t 
    cubicSolve(const double& p, const double& q, const double& r, 
	       double& root1, double& root2, double& root3)
    {
      static const double maxSqrt 
	= std::sqrt(std::numeric_limits<double>::max());
      static const double maxCubeRoot 
	= std::pow(std::numeric_limits<double>::max(), 1.0/3.0);

      if (r == 0)
	{
	  //no constant term, so divide by x and the result is a
	  //quadratic, but we must include the trivial x = 0 root
	  if (quadSolve(q,p,1.0, root1, root2))
	    {
	      root3 = 0;
	      if (root1 < root2) std::swap(root1,root2);
	      if (root2 < 0) 
		{
		  std::swap(root2,root3);
		  if (root1 < 0) std::swap(root1,root2);
		}
	      return 3;
	    }
	  else
	    {
	      root1 = 0;
	      return 1;
	    }
	}

      if ((p == 0) && (q == 0))
	{
	  //Special case
	  //Equation is x^3 == -r
	  root1 = std::pow(-r, 1.0 / 3.0);
	  return 1;
	}

      if ((p > maxSqrt) || (p < -maxSqrt))
	{
	  //Equation limits to x^3 + p * x^2 == 0
	  root1 = -p;
	  return 1;
	}

      if (q > maxSqrt)
	{
	  //Special case, if q is large and the root is -r/q,
	  //The x^3 term is negligble, and all other terms cancel.
	  root1 = -r / q;
	  return 1;
	}

      if (q < -maxSqrt)
	{
	  //Special case, equation is x^3 + q x == 0
	  root1 = -std::sqrt(-q);
	  return 1;
	}

      if ((r > maxSqrt) || (r < -maxSqrt))
	{
	  //Another special case
	  //Equation is x^3 == -r
	  root1 = std::pow(-r, 1.0 / 3.0);
	  return 1;
	}    

      //Need to insert the special cases for overflows here as well
      double u = q - p * p / 3.0;
      double v = r - p * q / 3.0 + 2.0 * p * p * p / 27.0;

      double j = 4.0 * (u / 3.0) * (u / 3.0) * (u / 3.0) + v * v;
  
      if (j >= 0) //Only one root
	{
	  double w = std::sqrt(j);
	  if (v < 0)
	    root1 = std::pow(0.5*(w-v), 1.0/3.0) - (u / 3.0) * std::pow(2.0 / (w-v), 1.0/3.0) - p / 3.0;
	  else
	    root1 = (u / 3.0) * std::pow(2.0 / (w+v), 1.0/3.0) - std::pow(0.5*(w+v), 1.0/3.0) - p / 3.0;
      
	  //Special cases for overflows
	  if (std::abs(p) > 27 * maxCubeRoot) 
	    root1 = - p;
	  if (std::abs(v) > maxSqrt) 
	    root1 = - std::pow(v, 1.0 / 3.0);
	  if (std::abs(u) > 0.75 * maxCubeRoot) 
	    root1 = std::pow(4, 1.0 / 3.0) * u / 3.0;
      
	  return 1;
	}
  
      double s = std::sqrt(-u / 3.0);
      double t = - v / (2.0 * s * s * s);
      double k = std::acos(t) / 3.0;
  
      root1 = 2 * s * std::cos(k) - p / 3.0;
      root2 = s * (-std::cos(k) + std::sqrt(3) * std::sin(k)) - p / 3.0;
      root3 = s * (-std::cos(k) - std::sqrt(3) * std::sin(k)) - p / 3.0;
      return 3;
    }
  }
}

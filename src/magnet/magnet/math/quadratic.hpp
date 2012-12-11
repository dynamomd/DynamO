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
#include <cmath>
#include <complex>
#include <magnet/exception.hpp>

namespace magnet {
  namespace math {
    /*! \brief An exception thrown when the quadratic equation has no
        roots, or the roots are complex and real-only form is used. 
    */
    class NoQuadraticRoots : public std::exception {};

    /*! \brief Solves a quadratic equation of the form
        \f$a\,x^2+b\,x+c=0\f$ and returns the (possibly complex) roots.
	
	\throw NoQuadraticRoots If \f$a=0\f$ and \f$b=0\f$ as this equation is
	not a function of \f$x\f$.
	
	\sa quadraticEquation

	\return The roots of the quadratic.
     */
    inline std::pair<std::complex<double>,std::complex<double> > 
    quadraticEquationComplex(const double a, const double b, const double c)
    {
      if (a == 0)
	{
	  if (b == 0) throw NoQuadraticRoots();
	  double root = - c / b;
	  return std::make_pair(std::complex<double>(root), 
				std::complex<double>(root));
	}
      
      double delta=(b * b - 4 * a * c);
      double inv_2a = 1 / (2 * a);
      double root = std::sqrt(std::abs(delta));
      double real = -b * inv_2a;
      double imag = root * inv_2a;

      if (delta >= 0)
	return std::make_pair(std::complex<double>(real - imag), 
			      std::complex<double>(real + imag));
      else
	return std::make_pair(std::complex<double>(real, -imag), 
			      std::complex<double>(real, imag));
    }

    /*! \brief Solves a quadratic equation of the form
        \f$a\,x^2+b\,x+c=0\f$ for the real roots.

	This implementation avoids a catastrophic cancellation of
	errors. See the following link for more details:
	http://en.wikipedia.org/wiki/Quadratic_equation#Floating_point_implementation
	
	It also handles the case when the polynomial being a linear
	function (\f$a=0\f$).

	\throw NoQuadraticRoots If \f$a=0\f$ and \f$b=0\f$ (this equation is
	not a function of \f$x\f$) or if the roots are complex.
	
	\sa quadraticEquationComplex

	\return The roots of the quadratic.
     */
    inline std::pair<double, double>
    quadraticEquation(const double a, const double b, const double c)
    {
      if (a == 0)
	{
	  if (b == 0) throw NoQuadraticRoots();
	  double root = - c / b;
	  return std::make_pair(root, root);
	}
      
      double discriminant = b * b - 4 * a * c;
      if (discriminant < 0) throw NoQuadraticRoots();
      double arg = std::sqrt(discriminant);
      double q = -0.5 * ( b + ((b < 0) ? -arg : arg));
      
      return std::make_pair(q / a, c / q);
    }

    inline bool quadSolve(const double& C, const double& B, const double& A, 
			  double& root1, double& root2)
    {
      // Contingency: if A = 0, not a quadratic = linear
      if(A == 0)
	{
	  //If B is zero then we have a NaN
	  if(B == 0) return false;
      
	  root1 = -1.0 * C / B;
	  root2 = root1;
	}

      double discriminant = (B * B) - (4 * A * C);
      
      //Cannot do imaginary numbers, yet
      if (discriminant < 0) return false;
      
      //This avoids a cancellation of errors. See
      //http://en.wikipedia.org/wiki/Quadratic_equation#Floating_point_implementation
      double t = -0.5 * ( B + ((B < 0) ? -1 : 1) * std::sqrt(discriminant));
      
      root1 = t / A;
      root2 = C / t;

      return true;
    }
  }
}

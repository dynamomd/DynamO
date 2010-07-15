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
#include "../base/constants.hpp"
#include "../base/is_exception.hpp"
#include <boost/array.hpp>
#include <algorithm>

struct magSortClass 
{ inline bool operator()(Iflt i, Iflt j) { return std::abs(i) < std::abs(j);} };  

inline 
long int rintfunc(const double& x)
{ return lrint(x); }

inline 
long int rintfunc(const float& x)
{ return lrintf(x); }


template<int X, int Y>
struct ctime_pow {
  static const int result = X * ctime_pow<X, Y-1>::result;
};

template<int X>
struct ctime_pow<X,1> {
  static const int result = X;
};

typedef enum
  {
    ROOT_SMALLEST_EITHER   =   1,
    ROOT_SMALLEST_POSITIVE =   2,
    ROOT_SMALLEST_NEGATIVE =   4,
    ROOT_LARGEST_EITHER    =   8,
    ROOT_LARGEST_POSITIVE  =  16,
    ROOT_LARGEST_NEGATIVE  =  32
  } rootTypeEnum;

//Solve a quadratic of the form A x^2 + B x + C == 0
inline bool quadSolve(const Iflt& C, const Iflt& B, const Iflt& A, Iflt& root1, Iflt& root2)
{
  // Contingency: if A = 0, not a quadratic = linear
  if(A == 0)
    {
      //If B is zero then we have a NaN
      if(B == 0) return false;
      
      root1 = -1.0 * C / B;
      root2 = root1;
    }
  else
  {
    Iflt discriminant = (B * B) - (4 * A * C);

    //Cannot do imaginary numbers, yet
    if (discriminant < 0) return false;
    
    //This avoids a cancellation of errors. See
    //http://en.wikipedia.org/wiki/Quadratic_equation#Floating_point_implementation
    Iflt t((B < 0)
	   ? -0.5 * (B-std::sqrt(discriminant))
	   : -0.5 * (B+std::sqrt(discriminant)));
    
    root1 = t / A;
    root2 = C / t;
  }

  return true;
}

template<rootTypeEnum rootType>
inline bool quadSolve(const Iflt& C, const Iflt& B, const Iflt& A, Iflt& ans)
{
  Iflt root1(0), root2(0);

  if (!quadSolve(C,B,A,root1,root2)) return false;

  switch (rootType)
    {
    case ROOT_SMALLEST_EITHER:
      ans = (fabs(root1) < fabs(root2)) ? root1 : root2;
      break;
    case ROOT_LARGEST_EITHER:
      ans = (fabs(root1) < fabs(root2)) ? root2 : root1;
    case ROOT_LARGEST_NEGATIVE:
      if (root1 < 0 && root2 < 0)
	ans = ((root1 < root2) ? root1 : root2);
      else if (root1 < 0)
	ans = root1;
      else if (root2 < 0)
	ans = root2;
      else
	return false;
      break;
    case ROOT_SMALLEST_NEGATIVE:
      if (root1 < 0 && root2 < 0)
	ans = ((root1 < root2) ? root2 : root1);
      else if (root1 < 0)
	ans = root1;
      else if (root2 < 0)
	ans = root2;
      else
	return false;
      break;
    case ROOT_LARGEST_POSITIVE:
      if (root1 > 0 && root2 > 0)
	ans = ((root1 > root2) ? root1 : root2);
      else if (root1 > 0)
	ans = root1;
      else if (root2 > 0)
	ans = root2;
      else
	return false;
      break;
    case ROOT_SMALLEST_POSITIVE:
      if (root1 > 0 && root2 > 0)
	ans = ((root1 > root2) ? root2 : root1);
      else if (root1 > 0)
	ans = root1;
      else if (root2 > 0)
	ans = root2;
      else
	return false;
      break;
    default:
      D_throw() << "Unknown root selected";
      break;
    }
  return true;
}

//Please read  http://linus.it.uts.edu.au/~don/pubs/solving.html
//For solving cubics like x^3 + p * x^2 + q * x + r == 0
//This function always returns root1 >= root2 >= root3!
inline size_t 
cubicSolve(const Iflt& p, const Iflt& q, const Iflt& r, Iflt& root1, Iflt& root2, Iflt& root3)
{
  static const Iflt maxSqrt = std::sqrt(std::numeric_limits<Iflt>::max());
  static const Iflt maxCubeRoot = std::pow(std::numeric_limits<Iflt>::max(), 1.0/3.0);


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
  Iflt u = q - p * p / 3.0;
  Iflt v = r - p * q / 3.0 + 2.0 * p * p * p / 27.0;

  Iflt j = 4.0 * (u / 3.0) * (u / 3.0) * (u / 3.0) + v * v;
  
  if (j > 0) //Only one root
    {
      Iflt w = std::sqrt(j);
      if (v < 0)
  	root1 = std::pow(0.5*(w-v), 1.0/3.0) - (u / 3.0) * std::pow(2.0 / (w-v), 1.0/3.0) - p / 3.0;
      else
  	root1 = (u / 3.0) * std::pow(2.0 / (w+v), 1.0/3.0) - std::pow(0.5*(w+v), 1.0/3.0) - p / 3.0;
      
      //Special cases for overflows
      if (std::abs(p) > 27 * maxCubeRoot) root1 = - p;
      if (std::abs(v) > maxSqrt) root1 = - std::pow(v, 1.0 / 3.0);
      if (std::abs(u) > 0.75 * maxCubeRoot) root1 = std::pow(4, 1.0 / 3.0) * u / 3.0;
      
      return 1;
    }
  
  Iflt s = std::sqrt(-u / 3.0);
  Iflt t = - v / (2.0 * s * s * s);
  Iflt k = std::acos(t) / 3.0;
  
  root1 = 2 * s * std::cos(k) - p / 3.0;
  root2 = s * (-std::cos(k) + std::sqrt(3) * std::sin(k)) - p / 3.0;
  root3 = s * (-std::cos(k) - std::sqrt(3) * std::sin(k)) - p / 3.0;
  return 3;
}

size_t 
neumarkQuarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d, 
		    Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4);

size_t
yacfraidQuarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d,
		     Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4);

//Solves quartics of the form x^4 + a x^3 + b x^2 + c x + d ==0
inline size_t quarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d, 
			   Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4)
{
  static const Iflt maxSqrt = std::sqrt(std::numeric_limits<Iflt>::max());

  if (std::abs(a) > maxSqrt)
    yacfraidQuarticSolve(a,b,c,d,root1,root2,root3,root4);

 if (d == 0)
    {//Solve a cubic with a trivial root of 0
      
      root1 = 0;
      return 1 + cubicSolve(a, b, c, root2, root3, root4);
    }
  
  if ((a == 0) && (c== 0))
    {//We have a biquadratic
      
      Iflt quadRoot1,quadRoot2;
      if (quadSolve(d,b,1, quadRoot1, quadRoot2))
	{
	  if (quadRoot1 < quadRoot2) std::swap(quadRoot1,quadRoot2);
	  
	  if (quadRoot1 < 0)
	    return 0;
	  
	  root1 = std::sqrt(quadRoot1);
	  root2 = -std::sqrt(quadRoot1);
	  
	  if (quadRoot2 < 0)
	    return 2;

	  root3 = std::sqrt(quadRoot2);
	  root4 = -std::sqrt(quadRoot2);
	  return 4;
	}
      else
	return 0;
      
    }


  //Now we have to resort to some dodgy formulae!
  return neumarkQuarticSolve(a,b,c,d,root1,root2,root3,root4);
  
}

Iflt quarticError(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d,
		  const boost::array<Iflt, 4>& roots, const size_t rootCount)
{
  boost::array<Iflt, 4> errors;

  for (size_t root = 0; root < rootCount; ++ root)
    {
      const Iflt value = (((roots[root]+a) * roots[root] + b) * roots[root] + c) * roots[root] + d;

      if (value == 0) { errors[root] = 0; continue; }

      const Iflt deriv = ((4 * roots[root] + 3 * a) * roots[root] + 2 * b) * roots[root] + c;
      
      if (deriv != 0) 
	errors[root] = std::abs(value / deriv);
      else
	{
	  const Iflt secDeriv = (12 * roots[root] + 6 * a) * roots[root] + 2 * b;
	  if (secDeriv != 0)
	    errors[root] = std::sqrt(std::abs(value / secDeriv));
	  else
	    {
	      const Iflt thirdDeriv = 24 * roots[root] + 6 * a;
	      if (thirdDeriv != 0)
		errors[root] = std::pow(std::abs(value / thirdDeriv), 1.0/3.0);
	      else
		errors[root] = std::sqrt(std::sqrt(std::abs(value)/24));
	    }
	}
    }

  return *std::min_element(errors.begin(), errors.begin()+rootCount);
}

//solve the quartic equation -
//x**4 + a*x**3 + b*x**2 + c*x + d = 0
size_t 
neumarkQuarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d, 
		    Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4)
{
  boost::array<Iflt, 4> rts;
  boost::array<Iflt, 3> worst3;
  int j,k, n4[4];
  double y,g,gg,h,hh,gdis,gdisrt,hdis,hdisrt,g1,g2,h1,h2;
  double bmy,gerr,herr,y4,bmysq;
  double v1[4],v2[4],v3[4];
  double hmax,gmax;
  double qrts[4][3];        /* quartic roots for each cubic root */

  if (d == 0.0)
    {
      root1 = 0.0;
      return cubicSolve(a,b,c,root2,root3,root4) + 1;
    }

  Iflt asq = a * a;

  Iflt d4 = d * 4.0;
  Iflt p =  -b * 2.0;
  Iflt q = b * b + a * c - d4;
  Iflt r = (c - a*b)*c + asq*d;
  size_t cubicRoots = cubicSolve(p,q,r,v3[0],v3[1],v3[2]);

  size_t nQuarticRoots[3]; 

  for (size_t j3 = 0; j3 < cubicRoots; ++j3)
    {
      y = v3[j3];
    
      bmy = b - y;
      y4 = y * 4.0;
      bmysq = bmy*bmy;
      gdis = asq - y4;
      hdis = bmysq - d4;

      if ((gdis < 0.0) || (hdis < 0.0))
	nQuarticRoots[j3] = 0;
      else
	{
	  g1 = a * 0.5;
	  h1 = bmy* 0.5;
	  gerr = asq + y4;
	  herr = hdis;
	  if (d > 0.0)
	    herr = bmysq + d4;
	  if ((y < 0.0) || (herr*gdis > gerr*hdis))
	    {
	      gdisrt = sqrt(gdis);
	      g2 = gdisrt*0.5;
	      if (gdisrt != 0.0)
		h2 = (a*h1 - c)/gdisrt;
	      else
		h2 = 0.0;
	    }
	  else
	    {
	      hdisrt = std::sqrt(hdis);
	      h2 = hdisrt*0.5;
	      if (hdisrt != 0.0)
		g2 = (a*h1 - c)/hdisrt;
	      else
		g2 = 0.0;
	    }
	  /*
	    note that in the following, the tests ensure non-zero
	    denominators -
	  */
	  h = h1 - h2;
	  hh = h1 + h2;
	  hmax = hh;
	  if (hmax < 0.0)
	      hmax =  -hmax;
	  if (hmax < h)
	      hmax = h;
	  if (hmax <  -h)
	      hmax =  -h;
	  if ((h1 > 0.0)&&(h2 > 0.0))
	      h = d/hh;
	  if ((h1 < 0.0)&&(h2 < 0.0))
	      h = d/hh;
	  if ((h1 > 0.0)&&(h2 < 0.0))
	      hh = d/h;
	  if ((h1 < 0.0)&&(h2 > 0.0))
	      hh = d/h;
	  if (h > hmax)
	      h = hmax;
	  if (h <  -hmax)
	      h =  -hmax;
	  if (hh > hmax)
	      hh = hmax;
	  if (hh < -hmax)
	      hh =  -hmax;

	  g = g1 - g2;
	  gg = g1 + g2;
	  gmax = gg;
	  if (gmax < 0.0)
	      gmax = -gmax;
	  if (gmax < g)
	      gmax = g;
	  if (gmax < -g)
	      gmax = -g;
	  if ((g1 > 0.0)&&(g2 > 0.0))
	      g = y/gg;
	  if ((g1 < 0.0)&&(g2 < 0.0))
	      g = y/gg;
	  if ((g1 > 0.0)&&(g2 < 0.0))
	      gg = y/g;
	  if ((g1 < 0.0)&&(g2 > 0.0))
	      gg = y/g;
	  if (g > gmax)
	      g = gmax;
	  if (g <  -gmax)
	      g = -gmax;
	  if (gg > gmax)
	      gg = gmax;
	  if (gg <  -gmax)
	      gg = -gmax;

	  size_t n1 = quadSolve(hh,gg, 1.0,v1[0],v1[1]);
	  size_t n2 = quadSolve(h,g,1.0,v2[0],v2[1]);
	  nQuarticRoots[j3] = 2*n1 + 2*n2;
	  qrts[0][j3] = v1[0];
	  qrts[1][j3] = v1[1];
	  qrts[2*n1+0][j3] = v2[0];
	  qrts[2*n1+1][j3] = v2[1];
	}

      for (j = 0; j < nQuarticRoots[j3]; ++j)
        rts[j] = qrts[j][j3];
      worst3[j3] = quarticError(a, b, c, d, rts,nQuarticRoots[j3]);

    } /* j3 loop */
size_t  j3 = 0;
 if (cubicRoots > 1)
  {
    if ((nQuarticRoots[1] > nQuarticRoots[j3]) ||
        ((worst3[1] < worst3[j3] ) && (nQuarticRoots[1] == nQuarticRoots[j3]))) j3 = 1;
    if ((nQuarticRoots[2] > nQuarticRoots[j3]) ||
        ((worst3[2] < worst3[j3] ) && (nQuarticRoots[2] == nQuarticRoots[j3]))) j3 = 2;
  }

  root1 = qrts[0][j3];
  root2 = qrts[1][j3];
  root3 = qrts[2][j3];
  root4 = qrts[3][j3];

  return nQuarticRoots[j3];
} /* neumark */
/****************************************************/

size_t 
descartesQuarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d, 
		      Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4)
{
  boost::array<Iflt, 4> rts;
  boost::array<Iflt, 3> worst3;
  double qrts[4][3];        /* quartic roots for each cubic root */

  int j, n4[4];
  double v1[4],v2[4],v3[4];
  double k,y;
  double dis;
  double p,q,r;
  double e0,e1,e2;
  double g,h;
  double asq;
  double ainv4;
  double e1invk;

  asq = a*a;
  e2 = b - asq * (3.0/8.0);
  e1 = c + a*(asq*0.125 - b*0.5);
  e0 = d + asq*(b*0.0625 - asq*(3.0/256.0)) - a*c*0.25;

  p = 2.0*e2;
  q = e2*e2 - 4.0*e0;
  r = -e1*e1;

  size_t n3 = cubicSolve(p,q,r,v3[0],v3[1],v3[2]);
  for (size_t j3 = 0; j3 < n3; ++j3)
  {
     y = v3[j3];
     if (y <= 0.0)
       n4[j3] = 0;
     else
     {
       k = std::sqrt(y);
       ainv4 = a*0.25;
       e1invk = e1/k;
       g = (y + e2 + e1invk)*0.5;
       h = (y + e2 - e1invk)*0.5 ;
       bool n1 = quadSolve( g, -k, 1.0, v1[0], v1[1]);
       bool n2 = quadSolve( h, k, 1.0, v2[0], v2[1]);
       qrts[0][j3] = v1[0] - ainv4;
       qrts[1][j3] = v1[1] - ainv4;
       qrts[n1*2][j3] = v2[0] - ainv4;
       qrts[n1*2+1][j3] = v2[1] - ainv4;
       n4[j3]= n1*2 + n2*2;
     } /* y>=0 */
donej3:
     for (j = 0; j < n4[j3]; ++j)
        rts[j] = qrts[j][j3];
     worst3[j3] = quarticError(a, b, c, d, rts, n4[j3]);
  } /* j3 loop */
done:
  size_t j3 = 0;
  if (n3 != 1)
  {
     if ((n4[1] > n4[j3]) ||
        ((worst3[1] < worst3[j3] ) && (n4[1] == n4[j3]))) j3 = 1;
     if ((n4[2] > n4[j3]) ||
        ((worst3[2] < worst3[j3] ) && (n4[2] == n4[j3]))) j3 = 2;
  }

  root1 = qrts[0][j3];
  root2 = qrts[1][j3];
  root3 = qrts[2][j3];
  root4 = qrts[3][j3];

  return (n4[j3]);
} /* descartes */
/****************************************************/


//solve the quartic equation -
//x**4 + a*x**3 + b*x**2 + c*x + d = 0
size_t 
yacfraidQuarticSolve(const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d,
		     Iflt& root1, Iflt& root2, Iflt& root3, Iflt& root4)
{
  int i,j;
  double y;
  double v1[4],v2[4],v3[4];
  double det0,det1,det2,det3;
  double det0rt,det1rt,det2rt,det3rt;
  double e,f,g,h,k;
  double fsq,gsq,hsq,invk;
  double P,Q,R,U;

  boost::array<Iflt, 4> rts;
  boost::array<Iflt, 3> worst3;
  double qrts[4][3];        /* quartic roots for each cubic root */

  if (d == 0.0)
    {
      root1 = 0.0;
      return cubicSolve(a,b,c,root2,root3,root4) + 1;
    }

  Iflt asq = a * a;
  Iflt acu = a * asq;
  Iflt b4 = b * 4.0;
  size_t n3 = 0;
  int n4[4];

  P = asq * b - b4 * b + 2.0 * a * c + 16.0 * d ;
  Q = asq * c - b4 * c + 8.0 * a * d;
  R = asq * d - c * c;
  U = acu - b4 * a + 8.0 * c;
  n4[0] = 0;


  asq = a*a;
  acu = a*asq;
  b4 = b*4.0;
  n3 = 0;

  P = asq*b - b4*b + 2.0*a*c + 16.0*d ;
  Q = asq*c - b4*c + 8.0*a*d;
  R = asq*d - c*c ;
  U = acu - b4*a + 8.0*c;
  n4[0] = 0;
  if (U == 0.0)
  {
     if (P == 0.0)
     {
        det0 = 3.0*asq - 8.0*b;
        if (det0 < 0.0)
           goto done;
        det0rt = sqrt(det0);
        qrts[0][0] = (-a + det0rt)*0.25;
        qrts[1][0] = qrts[0][0];
        qrts[2][0] = (-a - det0rt)*0.25;
        qrts[3][0] = qrts[2][0];
        n4[0] = 4;
        goto done;
     } /* P=0 */
     else
     {
        det1 = asq*asq - 8.0*asq*b + 16.0*b*b - 64.0*d;
        if (det1 < 0.0)
	  goto done;;
        n4[0] = 0;
        det1rt =  sqrt(det1);
        det2 = 3.0*asq - 8.0*b + 2.0*det1rt;
        if (det2 >= 0.0)
        {
           det2rt = sqrt(det2);
           qrts[0][0] = (-a + det2rt)*0.25;
           qrts[1][0] = (-a - det2rt)*0.25;
           n4[0] = 2;
        }
        det3 = 3.0*asq - 8.0*b - 2.0*det1rt;
        if (det3 >= 0.0)
        {
           det3rt = sqrt(det3);
           qrts[n4[0]++][0] = (-a + det3rt)*0.25;
           qrts[n4[0]++][0] = (-a - det3rt)*0.25;
        }
        goto done;
     } /* P<>0 */
  }

  n3 = cubicSolve(P/U,Q/U,R/U,v3[0], v3[1], v3[2]);
  for (size_t j3 = 0; j3 < n3; ++j3)
  {
     y = v3[j3];
     j = 0;
     k = a + 4.0*y;
     if (k == 0.0)
       goto donej3;
     invk = 1.0/k;
     e = (acu - 4.0*c - 2.0*a*b + (6.0*asq - 16.0*b)*y)*invk;
     fsq = (acu + 8.0*c - 4.0*a*b)*invk;
     if (fsq < 0.0)
       goto donej3;
     f = sqrt(fsq);
     gsq = 2.0*(e + f*k);
     hsq = 2.0*(e - f*k);
     if (gsq >= 0.0)
     {
        g = sqrt(gsq);
        qrts[j++][j3] = (-a - f - g)*0.25;
        qrts[j++][j3] = (-a - f + g)*0.25;
     }
     if (hsq >= 0.0)
     {
        h = sqrt(hsq);
        qrts[j++][j3] = (-a + f - h)*0.25;
        qrts[j++][j3] = (-a + f + h)*0.25;
     }
donej3:
     n4[j3] = j;
     for (j = 0; j < n4[j3]; ++j)
        rts[j] = qrts[j][j3];
     
     worst3[j3] = quarticError(a, b, c, d, rts, n4[j3]);
  } /* j3 loop */
done:
  size_t j3 = 0;
  if (n3 > 1)
  {
     if ((n4[1] > n4[j3]) ||
        ((worst3[1] < worst3[j3] ) && (n4[1] == n4[j3]))) j3 = 1;

     if ((n4[2] > n4[j3]) ||
        ((worst3[2] < worst3[j3] ) && (n4[2] == n4[j3]))) j3 = 2;
  }

  root1 = qrts[0][j3];
  root2 = qrts[1][j3];
  root3 = qrts[2][j3];
  root4 = qrts[3][j3];


  return (n4[j3]);
} /* yacfraid */
/*****************************************/

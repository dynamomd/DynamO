/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
namespace magnet {
  namespace math {
    class Spline : public std::vector<std::pair<double, double> >
    { 
      struct SplineData
      {
	double x;
	double a;
	double b;
	double c;
	double d;
      };

      std::vector<SplineData> _data;

    public:
      inline void addPoint(double x, double y)
      { push_back(std::pair<double, double>(x,y)); }

      void generate()
      {
	std::sort(begin(), end());

	const size_t e = size() - 1;

	ublas::matrix<double> A(size(), size());
	A(0,0) = 2 * h(0);
	A(1,0) = h(0);
	for (size_t i(1); i < e; ++i)
	  {
	    A(i-1,i) = h(i-1);
	    A(i,i) = 2 * (h(i-1) + h(i));
	    A(i+1,i) = h(i);
	  }
	A(e,e) = 2 * h(e - 1);
	A(e-1,e) = h(e - 1);

	ublas::vector<double> C(size());
	C(0) = (y(1) - y(0)) / h(0);
	C(e) = -(y(e) - y(e-1)) / h(e-1);
	for (size_t i(1); i < e; ++i)
	  C(i) = (y(i+1) - y(i)) / h(i)
	    - (y(i) - y(i-1)) / h(i-1);
	for (size_t i(0); i < size(); ++i)
	  C(i) *= 6;
	ublas::matrix<double> AInv(size(), size());
	InvertMatrix(A,AInv);
	
	_B = ublas::prod(C, AInv);

	_data.resize(size() - 1);
	for (size_t i(0); i < e; ++i)
	  {
	    _data[i].x = x(i);
	    _data[i].a = (_B(i+1) - _B(i)) / (6 * h(i));
	    _data[i].b = _B(i) / 2;
	    _data[i].c = (y(i+1) - y(i)) / h(i) - _B(i+1) * h(i) / 6 - _B(i) * h(i) / 3;
	    _data[i].d = y(i);
	  }
      }
      
      double operator()(double xval) const
      {
	if (xval <= x(0)) return y(0);
	if (xval >= x(size()-1)) return y(size()-1);

	for (std::vector<SplineData>::const_iterator iPtr = _data.begin();
	     iPtr != _data.end()-1; ++iPtr)
	    if ((xval >= iPtr->x) && (xval <= (iPtr+1)->x))
	      {
		double lx = xval - iPtr->x;
		
		return ((iPtr->a * lx + iPtr->b) * lx + iPtr->c) * lx + iPtr->d;
	      }

	double lx = xval - _data.back().x;
	return ((_data.back().a * lx + _data.back().b) * lx + _data.back().c) * lx + _data.back().d;
      }

    private:
      inline double x(size_t i) const { return operator[](i).first; }
      inline double y(size_t i) const { return operator[](i).second; }
      inline double h(size_t i) const { return x(i+1) - x(i); }

      template<class T>
      bool InvertMatrix(ublas::matrix<T> A,
			ublas::matrix<T>& inverse) 
      {
 	using namespace ublas;
	
 	// create a permutation matrix for the LU-factorization
 	permutation_matrix<std::size_t> pm(A.size1());
	
 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return false;
	
 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));
	
 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
	
 	return true;
      }

      ublas::vector<double> _B;
    };
  }
}

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
#include <algorithm>
#include <tr1/array>
#include <stddef.h>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

const size_t NDIM(3);

namespace magnet {
  namespace math {
    namespace {
      template<typename T> inline T SQR(T a) { return a * a; }
    };

    // list of operators that can act of Vectors and/or Matrices
    namespace ops {
      enum {
	NoOp,         // tag to indicate no operation
	RowOp,        // tag to indicate row operation
	ColumnOp,     // tag to indicate columns operation
	PlusOp,       // tag for addition operator    
	MinusOp,      // tag for subtraction operator
	TimesOp,      // tag for multiplication operator
	NegativeOp,   // tag for negation operator
	TransposeOp,  // tag to indicate transpose
	DyadicOp      // tag for dyadic operator
      };
    }

    // define vector expressions and matrix expressions 
    class Base;  
    template<class A=Base,int B=ops::NoOp,class C=Base> class VectorExpression;
    template<class A=Base,int B=ops::NoOp,class C=Base> class MatrixExpression;

    // 'typedef' the defaults as the Vector and Matrix types
    typedef VectorExpression<> Vector; 
    typedef MatrixExpression<> Matrix;

    // define default VectorExpression, which is the Vector type
    template <>
    class VectorExpression<>
    {
    public:
      double _data[3];// elements of the vector

      inline VectorExpression<>() {}             // default constructor

      // passive access to the elements through eval
      template<int I> inline double eval() const { return _data[I]; }

      inline VectorExpression (double a, double b, double c)
      { operator()(0) = a; operator()(1) = b; operator()(2) = c; }

      template<class T>
      inline VectorExpression(const std::tr1::array<T, 3>& vec)
      { for (size_t i(0); i < 3; ++i) operator()(i) = vec[i]; }

      /*! \brief Constuctor from a vector expression.
       */
      template<class A,int B,class C> 
      inline VectorExpression(const VectorExpression<A, B, C> &e)
      { for (size_t i(0); i < 3; ++i) operator()(i) = e(i); }

      // assign zero to all elements
      inline void zero() 
      { 
	for (size_t i(0); i < 3; ++i) 
	  operator()(i) = 0;
      }

      inline double nrm2() const 
      { 
	double sum(0);
	for (size_t i(0); i < 3; ++i) 
	  sum += operator()(i) * operator()(i);
	return sum; 
      }
      inline double nrm() const
      {
	double biggest = maxElement();
	double sqrbiggest = biggest * biggest;
	double invsqrbiggest = 1.0 / sqrbiggest;
	if (biggest!=0)
	  {
	    double sum(0);
	    for (size_t i(0); i < 3; ++i) 
	      sum += invsqrbiggest * operator()(i) * operator()(i);
	    return biggest * std::sqrt(sum);
	  }
	else
	  return 0;
      }

      inline double maxElement() const { return std::max(std::max(std::abs(_data[0]), std::abs(_data[1])), std::abs(_data[2])); }

      inline double& operator()(size_t i) { return _data[i]; }
      inline const double& operator()(size_t i) const { return _data[i]; }
      inline double& operator[](size_t i) { return operator()(i); }
      inline const double& operator[](size_t i) const { return operator()(i); }

      inline bool operator==(const VectorExpression<>& ovec) const
      {
	for (size_t i(0); i < 3; ++i)
	  if (operator()(i) != ovec(i))
	    return false;
	return true;
      }

      inline bool operator!=(const VectorExpression<>& ovec) const
      {	return !operator==(ovec); }


      /*! \brief Assignment of vector from a vector expression.
       */
      template<class A,int B,class C>
      inline VectorExpression<>& operator=(const VectorExpression<A,B,C> &e)
      { 
	/* 
	   Here, we take a copy of the expression result to stop
	   expressions with *this vector on the RHS from having a race
	   condition.
	 */
	double newvals[3] = {e.eval<0>(), e.eval<1>(), e.eval<2>()};

	for (size_t i(0); i < 3; ++i) operator()(i) = newvals[i];
	return *this; 
      }
  
      // addto-assignment from an expression
      template<class A,int B,class C> 
      inline VectorExpression<>& operator+=(const VectorExpression<A,B,C> &e)
      { return (*this = *this + e); }


      // subtract-from assignment from an expression
      template<class A,int B,class C> 
      inline VectorExpression<>& operator-=(const VectorExpression<A,B,C> &e)
      { return (*this = *this - e); }

      // multiply-by assignment from an expression
      inline VectorExpression<>& operator*=(const double d)
      { 
	for (size_t i(0); i < 3; ++i) operator()(i) *= d;
	return *this; 
      }
  
      // divide-by assignment from an expression
      inline VectorExpression<>& operator/=(const double d)
      { return (*this *= (1.0/d)); }
      
      std::string toString() const
      {
	std::ostringstream os;
	os << "<" << _data[0] << "," << _data[1] << "," << _data[2] << ">";
	return os.str();
      }
    };

    // definition of vector-matrix operations:

    // define body of const members function of VectorExpression

    // to access the elements of a vector expression
#define VECPARENTHESES operator()(const int _i) const	\
    {							\
      switch(_i){					\
      case 0: return eval<0>();				\
      case 1: return eval<1>();				\
      case 2: return eval<2>();				\
      default: return 0;				\
      }							\
    }

    // to get the norm of a vector
#define VECNRM nrm() const				\
    {							\
      double x=eval<0>();				\
      double y=eval<1>();				\
      double z=eval<2>();				\
      double biggest= std::max(std::max(std::abs(x), std::abs(y)), std::abs(z)); \
      if (biggest!=0)					\
	return biggest*(double)(sqrt(SQR(x/biggest)	\
				     +SQR(y/biggest)	\
				     +SQR(z/biggest))); \
      else						\
	return 0;					\
    }

    // to get the norm squared of a vector expression
#define VECNRM2 nrm2() const					\
    {								\
      return SQR(eval<0>()) + SQR(eval<1>()) + SQR(eval<2>());	\
    }

    // Expression template class for '+'  between two VectorExpressions
    //
    // A, C, D, and F are classes, B and E are Operators
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Vector
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Vector
    template<class A,int B,class C,class D,int E,class F>
    class VectorExpression<VectorExpression<A,B,C>,ops::PlusOp,VectorExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to
      template<int I> inline double eval() const
      { 
	return a->eval<I>() + b->eval<I>(); 
      }

      // constructor
      inline VectorExpression(const VectorExpression<A,B,C>& _a,
			      const VectorExpression<D,E,F>& _b)
	: a(&_a), b(&_b) {}

    private:
      // pointers to the sub-expressions
      const VectorExpression<A,B,C> * a;
      const VectorExpression<D,E,F> * b;
    };


    // Expression template class for '-'  between two VectorExpressions
    //
    // A, C, D, and F are classes, B and E are Operators
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Vector
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Vector
    template<class A,int B,class C,class D,int E,class F>
    class VectorExpression<VectorExpression<A,B,C>,ops::MinusOp,VectorExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to
      template <int I> inline double eval() const
      { 
	return a->eval<I>() - b->eval<I>(); 
      }

      // constructor
      inline VectorExpression(const VectorExpression<A,B,C>& _a,
			      const VectorExpression<D,E,F>& _b)
	: a(&_a), b(&_b) {}

    private:
      // pointers to the sub-expressions 
      const VectorExpression<A,B,C> * a;
      const VectorExpression<D,E,F> * b;
    };


    // Expression template class for '^' between two VectorExpressions,
    // i.e. their outer product
    //
    // A, C, D, and F are classes, B and E are Operators
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Vector
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Vector
    template<class A,int B,class C,class D,int E,class F>
    class VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, VectorExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES
 
      // define what this operation evaluates to
      template <const int I> inline double eval() const
      {
	// despite its looks, a switch statement is actually quite efficient
	switch (I) {
	case 0:
	  return a->eval<1>() * b->eval<2>() - a->eval<2>() * b->eval<1>();
	case 1:
	  return a->eval<2>() * b->eval<0>() - a->eval<0>() * b->eval<2>();
	case 2:
	  return a->eval<0>() * b->eval<1>() - a->eval<1>() * b->eval<0>();
	default: 
	  return 0; //should not happen
	}
      }

      // constructor
      inline VectorExpression(const VectorExpression<A,B,C>& _a,
			      const VectorExpression<D,E,F>& _b)
	: a(&_a), b(&_b) {}

    private:
      // pointers to the sub-expressions 
      const VectorExpression<A,B,C> * a;
      const VectorExpression<D,E,F> * b;
    };


    // Expression template class for '*' between a Vector and a double
    //
    // A and C are classes, B is an Operator
    // if B=ops::NoOp and A=C=Base, the VectorExpression is actually a Vector
    template<class A,int B,class C>
    class VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2  
      inline double VECNRM
      inline double VECPARENTHESES
 
      // define what this operation evaluates to
      template <int I> inline double eval() const
      { 
	return a->eval<I>() * b; 
      }

      // optimize a second multiplication with a double
      inline VectorExpression & operator* (double c)
      {
	b *= c;
	return *this;
      }

      // optimize a subsequent division by a double
      inline VectorExpression & operator/ (double c)
      {
	b /= c;
	return *this;
      }

      // constructor
      inline VectorExpression(const VectorExpression<A,B,C> & _a, 
			      double _b) 
	: a(&_a), b(_b) {}

    private:
      // pointer to the subexpression 
      const VectorExpression<A,B,C>* a;
      // the double to multiply with
      double b;

      // be-friend multiplication operators that optimize further double
      // multiplication
      template<class D,int E,class F> 
      friend VectorExpression < VectorExpression<D,E,F>, ops::TimesOp, double >&
      operator* (double b, VectorExpression < VectorExpression<D,E,F>, ops::TimesOp, double > &a);

      template<class D,int E,class F> 
      friend VectorExpression < VectorExpression<D,E,F>, ops::TimesOp, double >&
      operator/ (double b, VectorExpression < VectorExpression<D,E,F>, ops::TimesOp, double > &a);

    };


    // Expression template class for unary '-' acting on a VectorExpression
    //
    // A and C are classes, B is an Operator
    // if B=ops::NoOp and A=C=Base, the VectorExpression is actually a Vector
    template<class A,int B,class C>
    class VectorExpression<VectorExpression<A,B,C>, ops::NegativeOp, Base>
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to
      template <int I> inline double eval() const
      {
	return  - a->eval<I>();
      }

      // constructor
      inline VectorExpression(const VectorExpression<A,B,C> & _a) 
	: a(&_a) {}

    private:
      // store pointer to sub-expression
      const VectorExpression<A,B,C> * a;

    };
    //
    // definition of the operators
    //

    // Vector + Vector
    template<class A,int B,class C,class D,int E,class F> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::PlusOp, VectorExpression<D,E,F> >
    operator+ (const VectorExpression<A,B,C> & a, 
	       const VectorExpression<D,E,F> & b)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::PlusOp, VectorExpression<D,E,F> >(a,b);
    }

    // Vector - Vector
    template<class A,int B,class C,class D,int E,class F> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::MinusOp, VectorExpression<D,E,F> >
    operator- (const VectorExpression<A,B,C> & a, 
	       const VectorExpression<D,E,F> & b)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::MinusOp, VectorExpression<D,E,F> >(a, b);
    }

    // Vector ^ Vector
    template<class A,int B,class C,class D,int E,class F> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, VectorExpression< D,E,F> >
    operator^ (const VectorExpression<A,B,C> & a, 
	       const VectorExpression<D,E,F> & b)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, VectorExpression< D,E,F> > (a, b);
    }

    // double * Vector
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>
    operator* (double b, 
	       const VectorExpression<A,B,C> & a)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>(a, b);
    }

    // Vector * double
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>
    operator* (const VectorExpression<A,B,C> & a, 
	       double b)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>(a, b);
    }

    // Vector / double
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>
    operator/ (const VectorExpression<A,B,C> & a, 
	       double b)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>(a, 1.0/b);
    }

    // Vector * double * double
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>&
    operator*(double b,VectorExpression<VectorExpression<A,B,C>,ops::TimesOp,double>&a)
    { 
      a.b *= b; 
      return a; 
    }

    // Vector * double / double
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::TimesOp, double>&
    operator/(double b,VectorExpression<VectorExpression<A,B,C>,ops::TimesOp,double>&a)
    { 
      a.b /= b; 
      return a; 
    }

    // - Vector
    template<class A,int B,class C> 
    inline VectorExpression<VectorExpression<A,B,C>, ops::NegativeOp, Base>
    operator-(const VectorExpression<A,B,C> & a)
    { 
      return VectorExpression<VectorExpression<A,B,C>, ops::NegativeOp, Base> (a);
    }

    // Vector | Vector
    template<class A,int B,class C,class D,int E,class F>
    inline double 
    operator|(const VectorExpression<A,B,C> & a, 
	      const VectorExpression<D,E,F> & b)
    {
      return 
	a.template eval<0>()*b.template eval<0>() 
	+ a.template eval<1>()*b.template eval<1>() 
	+ a.template eval<2>()*b.template eval<2>();
    }

    // Vector * Vector
    template<class A,int B,class C,class D,int E,class F>
    inline double 
    operator*(const VectorExpression<A,B,C> & a, const VectorExpression<D,E,F> & b)
    {
      return 
	a.template eval<0>()*b.template eval<0>() 
	+ a.template eval<1>()*b.template eval<1>() 
	+ a.template eval<2>()*b.template eval<2>();
    }

    // vectors
    template<class A, int B, class C>
    inline magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
					      const VectorExpression<A,B,C> & t )
    {
      char name[2] = "x";
  
      for (size_t iDim = 0; iDim < NDIM; iDim++)
	{
	  name[0]= 'x'+iDim; //Write the dimension
	  XML << magnet::xml::attr(name) << t(iDim);
	}
  
      return XML;
    }

    inline
    VectorExpression<>& 
    operator<<(VectorExpression<>& data, const magnet::xml::Node& XML)
    {
      for (size_t iDim = 0; iDim < NDIM; iDim++) 
	{
	  char name[2] = "x";
	  name[0] = 'x' + iDim; //Write the name
	  if (!XML.getAttribute(name))
	    name[0] = '0'+iDim;
      
	  try {
	    data[iDim] = XML.getAttribute(name).as<double>();
	  }
	  catch (boost::bad_lexical_cast &)
	    {
	      M_throw() << "Failed a lexical cast in CVector";
	    }
	}

      return data;
    }
  }
}

namespace coil { typedef ::magnet::math::Vector Vector; }
namespace dynamo { typedef ::magnet::math::Vector Vector; }

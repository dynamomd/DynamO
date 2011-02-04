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

//
// vecmat3.h
//
// Fast three dimensional "Vector" and "Matrix" classes using templates
//
// Ramses van Zon
// May 7, 2007
//
//  - Uses a predefined "DOUBLE", which if not defined, is set to "double"
//  - See 'docvecmat3.tex for further documentation on how to use these classes
//
//
// On the implementations of the template expressions:
//
// To have e.g. 'c = a+b' be unrolled to 
// 'c.x = a.x+b.x; c.y = a.y+b.y; c.z = a.z+b.z;',
// an expression template class is defined with the following properties:
//
// - It contains pointers to the vectors 'a' AND 'b'.
//
// - The constructor of this class sets up these pointers
//
// - It has an integer template parameter indicating the type of operation.
//   The possible operations are listed in the unnamed 'enum' below.
//
// - It contains a function 'eval<I>()', where I is an integer parameter, 
//   0, 1 or 2, which returns the value of the Ith element of the vector 
//   expression. [For matrix expressions this is replaced by 'eval<I,J>()'.]
//
// - The 'operator+' is overloaded to return an instance of this class
//
// - The 'operator=' is overloaded to call the 'eval<I>()' function.
//
// - Two class template-parameters are needed to specify the operant types.
//   This way one can distinguish different operations, e.g. vector-vector 
//   multiplication from vector-double multiplication. 
//
// - There are separate expression templates Vector-valued and Matrix-valued 
//   expressions, such that one can distinguish e.g. Matrix-Vector from 
//   Matrix-Matrix multiplication.
//
// The general expression templates are thus:
//   template<class A, int B, class C>
//   VectorExpression<A,B,C>; 
// and
//   template<class A, int B, class C>
//   MatrixExpression<A,B,C>;
// 
// The general template classes are not defined, only special instances for 
// allowed operations B. 
//
// The default values for the templates arguments are 'Base' for 'A'
// and 'C', and 'NoOp' for 'B' The templates with these default
// values, denoted as 'VectorExpression<>' and 'MatrixExpression<>'
// serve as the actual Vector and Matrix class, and are 'typedef'ed as
// such.  The specializations of these default templates therefore
// contain the actual elements and much of the basic functionality of
// vectors and matrices.
//
// Note that the classes A and C can themselved be expression classes,
// and so nested expressions such as (a+b)*(d+2*c) are perfectly
// possible.  This recursiveness does make the notation somethat
// involved.
//
// For further information on the technique of expression templates in c++, see
//  - ubiety.uwaterloo.ca/~tveldhui/papers/Expression-Templates/exprtmpl.html
//  - www.oonumerics.org/blitz
//  - tvmet.sourceforge.net
 
#include <cmath>
#include <algorithm>

namespace {
  template<typename T> inline T SQR(T a) { return a * a; }
};

// list of operators that can act of Vectors and/or Matrices
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

// define vector expressions and matrix expressions 
class Base;  
template<class A=Base,int B=NoOp,class C=Base> class VectorExpression;
template<class A=Base,int B=NoOp,class C=Base> class MatrixExpression;

// 'typedef' the defaults as the Vector and Matrix types
typedef VectorExpression<> Vector; 
typedef MatrixExpression<> Matrix;

// define default VectorExpression, which is the Vector type
template <>
class VectorExpression<>
{
 public:
  double x;                                  // elements of the vector
  double y;
  double z;

  inline VectorExpression<>() {}             // default constructor

  inline VectorExpression (double a, double b, double c)
    : x(a), y(b), z(c) {}

  // constuctor from a vector expression
  template<class A,int B,class C> 
  inline VectorExpression(const VectorExpression<A, B, C> &e):
    x(e.eval<0>()),
    y(e.eval<1>()),
    z(e.eval<2>())
  {}

  // assign zero to all elements
  inline void zero() { x = y = z = 0;}

  inline double nrm2() const { return x*x+y*y+z*z; }
  inline double nrm() const
  {
    double biggest=std::max(std::max(std::abs(x), std::abs(y)), std::abs(z));
    if (biggest!=0)
      return biggest*(double)(std::sqrt(SQR(x/biggest)
					+SQR(y/biggest)
					+SQR(z/biggest)));
    else
      return 0;
  }

  inline double& operator()(int i) { return *(&x+i); }// return element i
  inline const double& operator()(int i) const { return *(&x+i); }// return element i

  inline double& operator[](int i) { return *(&x+i); }// return element i
  inline const double& operator[](int i) const { return *(&x+i); };

  // assignment from an expression
  template<class A,int B,class C> 
  inline VectorExpression<>& operator=(const VectorExpression<A,B,C> &e)
  { x = e.eval<0>(); y = e.eval<1>(); z = e.eval<2>(); return *this; }
  
  // addto-assignment from an expression
  template<class A,int B,class C> 
  inline VectorExpression<>& operator+=(const VectorExpression<A,B,C> &e)
  { x += e.eval<0>(); y += e.eval<1>(); z += e.eval<2>(); return *this; }


  // subtract-from assignment from an expression
  template<class A,int B,class C>
  inline VectorExpression<>& operator-=(const VectorExpression<A,B,C> &e)
  { x -= e.eval<0>(); y -= e.eval<1>(); z -= e.eval<2>(); return *this; }

  // multiply-by assignment from an expression
  inline VectorExpression<>& operator*=(const double d)
  { x *= d; y *= d; z *= d; return *this; }
  
  // divide-by assignment from an expression
  inline VectorExpression<>& operator/=(const double d)
  { x /= d; y /= d; z /= d; return *this; }

  // passive access to the elements through eval
  template<int I> inline double eval() const;// template evaluation
};

template <> inline double VectorExpression<>::eval<0>() const { return x; }
template <> inline double VectorExpression<>::eval<1>() const { return y; }
template <> inline double VectorExpression<>::eval<2>() const { return z; }

// definition of vector-matrix operations:

// operator+ takes two vector expressions and returns an appropriate
// vector expression
template<class A,int B,class C,class D,int E,class F> 
inline VectorExpression<VectorExpression<A,B,C>,PlusOp,VectorExpression<D,E,F> >
operator+(const VectorExpression<A,B,C> & a, 
	  const VectorExpression<D,E,F> & b);

// define operator- as takes two vector expressions and returns an
// appropriate vector expression
template<class A,int B,class C,class D,int E,class F> 
inline VectorExpression<VectorExpression<A,B,C>,MinusOp,VectorExpression<D,E,F> > 
operator-(const VectorExpression<A,B,C> & a, 
	  const VectorExpression<D,E,F> & b);

// operator^ as takes two vector expressions and returns an
// appropriate vector expression
template<class A,int B,class C,class D,int E,class F> 
inline VectorExpression<VectorExpression<A,B,C>,TimesOp,VectorExpression<D,E,F> > 
operator^(const VectorExpression<A,B,C> & a, 
	  const VectorExpression<D,E,F> & b);

// operator* takes a scalar and a vector expression, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator* (double b, 
	   const VectorExpression<A,B,C> & a);

// operator* takes a vector expression and a scalar, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator* (const VectorExpression<A,B,C> & a, 
	   double b);

// operator/ takes a vector expression and a scalar, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator/ (const VectorExpression<A,B,C> & a, 
	   double b);

// unary operator- takes a vector expression, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, NegativeOp, Base>
operator- (const VectorExpression<A,B,C> & a);

// define dot product as operator| as takes two vector expressions and
// returns a double
template<class A,int B,class C,class D,int E,class F>
inline double 
operator|(const VectorExpression<A,B,C> &a, 
	  const VectorExpression<D,B,C> &b);

// also define dot product as operator* as takes two vector
// expressions and returns a double
template<class A,int B,class C,class D,int E,class F>
inline double 
operator* (const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,B,C> & b);

// distance between two vectors
template<class A,int B,class C,class D,int E,class F>
inline double dist(const VectorExpression<A,B,C> &a, 
		   const  VectorExpression<D,E,F> &b); 

// distance squared between two vectors
template<class A,int B,class C,class D,int E,class F>
inline double dist2(const VectorExpression<A,B,C> &a, 
		    const  VectorExpression<D,E,F> &b);

// distance between two vectors with the first one shifted by s
template<class A,int B,class C,class D,int E,class F,class G,int H,class I>
inline double distwithshift(const VectorExpression<A,B,C> & a, 
			    const VectorExpression<D,E,F> & b,
			    const VectorExpression<G,H,I> & s);

//
// INLINE IMPLEMENTATION
//

//
// member functions of VectorExpression and MatrixExpression
//


// define body of const members function of VectorExpression

// to access the elements of a vector expression
#define VECPARENTHESES operator()(const int _i) const \
{ \
  switch(_i){ \
    case 0: return eval<0>(); \
    case 1: return eval<1>(); \
    case 2: return eval<2>(); \
    default: return 0; \
  } \
}

// to get the norm of a vector
#define VECNRM nrm() const\
{\
  double x=eval<0>(); \
  double y=eval<1>(); \
  double z=eval<2>(); \
  double biggest=x<0?-x:x; \
  double absval=y<0?-y:y; \
  if (absval>biggest) biggest=absval; \
  absval=z<0?-z:z; \
  if (absval>biggest) biggest=absval; \
  if (biggest!=0) \
    return biggest*(double)(sqrt(SQR(x/biggest) \
			       +SQR(y/biggest) \
			       +SQR(z/biggest))); \
  else \
    return 0; \
}

// to get the norm squared of a vector expression
#define VECNRM2 nrm2() const \
{\
   return SQR(eval<0>()) + SQR(eval<1>()) + SQR(eval<2>()); \
}

// Expression template class for '+'  between two VectorExpressions
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Vector
// if E=NoOp and D=F=Base, the second sub-expression is actually a Vector
template<class A,int B,class C,class D,int E,class F>
class VectorExpression<VectorExpression<A,B,C>,PlusOp,VectorExpression<D,E,F> >
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
// if B=NoOp and A=C=Base, the first sub-expression is actually a Vector
// if E=NoOp and D=F=Base, the second sub-expression is actually a Vector
template<class A,int B,class C,class D,int E,class F>
class VectorExpression<VectorExpression<A,B,C>,MinusOp,VectorExpression<D,E,F> >
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
// if B=NoOp and A=C=Base, the first sub-expression is actually a Vector
// if E=NoOp and D=F=Base, the second sub-expression is actually a Vector
template<class A,int B,class C,class D,int E,class F>
class VectorExpression<VectorExpression<A,B,C>, TimesOp, VectorExpression<D,E,F> >
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
// if B=NoOp and A=C=Base, the VectorExpression is actually a Vector
template<class A,int B,class C>
class VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
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
  friend VectorExpression < VectorExpression<D,E,F>, TimesOp, double >&
  operator* (double b, VectorExpression < VectorExpression<D,E,F>, TimesOp, double > &a);

  template<class D,int E,class F> 
  friend VectorExpression < VectorExpression<D,E,F>, TimesOp, double >&
  operator/ (double b, VectorExpression < VectorExpression<D,E,F>, TimesOp, double > &a);

};


// Expression template class for unary '-' acting on a VectorExpression
//
// A and C are classes, B is an Operator
// if B=NoOp and A=C=Base, the VectorExpression is actually a Vector
template<class A,int B,class C>
class VectorExpression<VectorExpression<A,B,C>, NegativeOp, Base>
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
inline VectorExpression<VectorExpression<A,B,C>, PlusOp, VectorExpression<D,E,F> >
operator+ (const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,E,F> & b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, PlusOp, VectorExpression<D,E,F> >(a,b);
}

// Vector - Vector
template<class A,int B,class C,class D,int E,class F> 
inline VectorExpression<VectorExpression<A,B,C>, MinusOp, VectorExpression<D,E,F> >
operator- (const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,E,F> & b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, MinusOp, VectorExpression<D,E,F> >(a, b);
}

// Vector ^ Vector
template<class A,int B,class C,class D,int E,class F> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, VectorExpression< D,E,F> >
operator^ (const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,E,F> & b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, VectorExpression< D,E,F> > (a, b);
}

// double * Vector
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator* (double b, 
	   const VectorExpression<A,B,C> & a)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, double>(a, b);
}

// Vector * double
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator* (const VectorExpression<A,B,C> & a, 
	   double b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, double>(a, b);
}

// Vector / double
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>
operator/ (const VectorExpression<A,B,C> & a, 
	   double b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, double>(a, 1.0/b);
}

// Vector * double * double
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>&
operator*(double b,VectorExpression<VectorExpression<A,B,C>,TimesOp,double>&a)
{ 
  a.b *= b; 
  return a; 
}

// Vector * double / double
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, double>&
operator/(double b,VectorExpression<VectorExpression<A,B,C>,TimesOp,double>&a)
{ 
  a.b /= b; 
  return a; 
}

// - Vector
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, NegativeOp, Base>
operator-(const VectorExpression<A,B,C> & a)
{ 
  return VectorExpression<VectorExpression<A,B,C>, NegativeOp, Base> (a);
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

//
// Non-member functions
//

// return |a-b|
template<class A,int B,class C,class D,int E,class F>
inline double dist(const VectorExpression<A,B,C> & a, 
		   const VectorExpression<D,E,F> & b)
{
  double x = fabs(a.template eval<0>()-b.template eval<0>());
  double y = fabs(a.template eval<1>()-b.template eval<1>());
  double z = fabs(a.template eval<2>()-b.template eval<2>());
  double biggest=x;
  if (y>biggest) biggest=y;
  if (z>biggest) biggest=z;   
  if (biggest!=0) {
    x /= biggest;
    y /= biggest;
    z /= biggest;
    return biggest*sqrt(x*x+y*y+z*z);
  }
  else 
    return 0;
}

// return (a-b)|(a-b)
template<class A,int B,class C,class D,int E,class F>
inline double dist2(const VectorExpression<A,B,C> & a, 
		    const VectorExpression<D,E,F> & b)
{
  double d = a.template eval<0>() - b.template eval<0>();
  double c = d * d;
  d = a.template eval<1>() - b.template eval<1>();
  c += d * d;
  d = a.template eval<2>() - b.template eval<2>();
  return c + d * d;
}

// return |a+s-b|
template<class A,int B,class C,class D,int E,class F,class G,int H,class I>
inline double distwithshift(const VectorExpression<A,B,C> &a, 
			    const VectorExpression<D,E,F> & b, 
			    const VectorExpression<G,H,I>& s)
{
 double x=fabs(s.template eval<0>()+a.template eval<0>()-b.template eval<0>());
 double y=fabs(s.template eval<1>()+a.template eval<1>()-b.template eval<1>());
 double z=fabs(s.template eval<2>()+a.template eval<2>()-b.template eval<2>());
 double biggest=x;
 if (y>biggest) biggest=y;
 if (z>biggest) biggest=z;   
 if (biggest!=0) {
   x /= biggest;
   y /= biggest;
   z /= biggest;
   return biggest*sqrt(x*x+y*y+z*z);
 }
 else 
   return 0;
}

//#undef VECNRM
//#undef VECNRM2
//#undef VECPARENTHESIS

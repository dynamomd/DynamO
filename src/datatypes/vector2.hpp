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
 
#ifndef _VECMAT3_
#define _VECMAT3_

#include <cstdlib>
#include <cmath>

// type of elements is set by #defining the macro DOUBLE (default = double)
#ifndef DOUBLE
#define DOUBLE double
#endif

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

// forward declaration needed for the definitions of
// VectorExpression<> and MatrixExpression<>
class CommaOp;

// define default VectorExpression, which is the Vector type
template <>
class VectorExpression<>
{
 public:
  DOUBLE x;                                  // elements of the vector
  DOUBLE y;
  DOUBLE z;

  inline VectorExpression<>() {}             // default constructor
  inline explicit VectorExpression<>         // from three numbers
    (DOUBLE a,DOUBLE b,DOUBLE c); 
    
  template<class A,int B,class C>            // construct from VectorExpression
  inline VectorExpression<>
  	(const VectorExpression<A,B,C>&e);

  inline void zero();                        // make zero

  inline DOUBLE nrm() const;                 // norm
  inline DOUBLE nrm2() const;                // square norm

  inline DOUBLE& operator()(int i);          // return element i
  inline const DOUBLE& operator()(int i) const;

  inline DOUBLE& operator[](int i);          // return element i
  inline const DOUBLE& operator[](int i) const;

  inline VectorExpression<>& operator*=(const DOUBLE d); // multiply by number
  inline VectorExpression<>& operator/=(const DOUBLE d); // divide by number

  inline CommaOp operator=(DOUBLE c);        // start comma expression

  template<class A,int B,class C> inline 
  VectorExpression<>& 
  operator=(const VectorExpression<A,B,C> & e); //assign VectorExpression
  template<class A,int B,class C> inline 
  VectorExpression<>& 
  operator+=(const VectorExpression<A,B,C> & e);//add VectorExpression
  template<class A,int B,class C> inline 
  VectorExpression<>& 
  operator-=(const VectorExpression<A,B,C> & e);//subtract VectorExpression

  template<int I> inline DOUBLE eval() const;// template evaluation
};


// define default MatrixExpression, which is the Matrix type
template <>
class MatrixExpression<>
{
 public:
  DOUBLE xx, xy, xz,                         // elements of the matrix
         yx, yy, yz,
         zx, zy, zz;

  inline MatrixExpression<>() {}             // default constructor
  inline explicit MatrixExpression<>         // from 9 numbers
    (DOUBLE a,     DOUBLE b = 0, DOUBLE c = 0,
     DOUBLE d = 0, DOUBLE e = 0, DOUBLE f = 0,
     DOUBLE g = 0, DOUBLE h = 0, DOUBLE i = 0);

  template<class A,int B,class C>            // construct from MatrixExpression
  inline MatrixExpression<>
    (const MatrixExpression<A,B,C>&e);

  inline void zero();                        // make zero
  inline void one();                         // make identity 
  inline void reorthogonalize();             // Gramm-Schmidt
  
  inline DOUBLE nrm2() const;                // norm squared
  inline DOUBLE nrm() const;                 // norm 
  inline DOUBLE tr() const;                  // trace 
  inline DOUBLE det() const;                 // determinant 

  static inline Matrix identity() { return Matrix(1,0,0,0,1,0,0,0,1); }

  inline DOUBLE& operator()(int i, int j);   // return element i,j
  inline const DOUBLE&operator()(int i,int j) const;

  inline VectorExpression<MatrixExpression<>,RowOp,Base> 
  row(int i);                                // return i-th row   
  
  inline VectorExpression<MatrixExpression<>,ColumnOp,Base> 
  column(int j);                             // return j-th column
  
  template<class A, int B, class C> inline 
  void setRow(int i,                         // assign vector to row i
	      const VectorExpression<A,B,C>&e);

  template<class A, int B, class C> inline 
  void setColumn(int j,                      // or to column j
		 const VectorExpression<A,B,C>&e); 

  inline MatrixExpression<>& operator*=(const DOUBLE d);  // multiply by number
  inline MatrixExpression<>& operator/=(const DOUBLE d);  // divide by number

  inline CommaOp operator=(DOUBLE c);        // start of comma expression

  template<class A,int B,class C> inline     // assign MatrixExpression
  MatrixExpression<>& 
  operator=(const MatrixExpression<A,B,C>& e);   

  template<class A,int B,class C> inline 
  MatrixExpression<>& 
  operator+=(const MatrixExpression<A,B,C>& e);//add a MatrixExpression

  template<class A,int B,class C> inline 
  MatrixExpression<>& 
  operator-=(const MatrixExpression<A,B,C> &e);//subtract MatrixExpression

  template<class A,int B,class C> inline 
  MatrixExpression<>& 
  operator*=(const MatrixExpression<A,B,C> &m);//multiply MatrixExpression

  template<int I,int J> inline DOUBLE eval() const; //evaluate components
};
 
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
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator* (DOUBLE b, 
	   const VectorExpression<A,B,C> & a);

// operator* takes a vector expression and a scalar, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator* (const VectorExpression<A,B,C> & a, 
	   DOUBLE b);

// operator/ takes a vector expression and a scalar, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator/ (const VectorExpression<A,B,C> & a, 
	   DOUBLE b);

// unary operator- takes a vector expression, and returns an
// appropriate vector expression
template<class A,int B,class C>
inline VectorExpression<VectorExpression<A,B,C>, NegativeOp, Base>
operator- (const VectorExpression<A,B,C> & a);

// define dot product as operator| as takes two vector expressions and
// returns a DOUBLE
template<class A,int B,class C,class D,int E,class F>
inline DOUBLE 
operator|(const VectorExpression<A,B,C> &a, 
	  const VectorExpression<D,B,C> &b);

// also define dot product as operator* as takes two vector
// expressions and returns a DOUBLE
template<class A,int B,class C,class D,int E,class F>
inline DOUBLE 
operator* (const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,B,C> & b);

// operator+ takes two matrix expressions and returns an appropriate
// matrix expression
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, PlusOp, MatrixExpression<D,E,F> >
operator+(const MatrixExpression<A,B,C> & a, 
	  const MatrixExpression<D,E,F> & b);

// operator- takes two matrix expressions and returns an appropriate
// matrix expression
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, MinusOp, MatrixExpression<D,E,F> >
operator- (const MatrixExpression<A,B,C> & a,
	   const MatrixExpression<D,E,F> & b);

// operator* takes two matrix expressions and returns an appropriate
// matrix expression
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, MatrixExpression<D,E,F> >
operator*(const MatrixExpression<A,B,C> & a, 
	  const MatrixExpression<D,E,F> & b);

// operator* takes a matrix expressions and a vector expression, and
// returns an appropriate vector expression
template<class A,int B,class C,class D,int E,class F>
inline VectorExpression < MatrixExpression<A,B,C>, TimesOp, VectorExpression<D,E,F> >
operator*(const MatrixExpression<A,B,C> & a, 
	  const VectorExpression<D,E,F> & b);

// operator* takes a scalar and a matrix expression, and returns an
// appropriate matrix expression
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator* (DOUBLE b, 
	   const MatrixExpression<A,B,C> & a);

// operator* takes a matrix expression and a scalar, and returns an
// appropriate matrix expression
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator* (const MatrixExpression<A,B,C> & a, 
	   DOUBLE b);

// operator/ takes a matrix expression and a scalar, and returns an
// appropriate matrix expression
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator/ (const MatrixExpression<A,B,C> & a, 
	   DOUBLE b);

// unary operator- takes a matrix expression and returns an
// appropriate matrix expression
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, NegativeOp, Base >
operator- (const MatrixExpression<A,B,C> & a);

// Transpose function, takes a matrix expression and returns an
// appropriate matrix expression
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TransposeOp, Base> 
Transpose(const MatrixExpression<A,B,C> & a);

// Dyadic function, takes two vector expressions and returns an
// appropriate matrix expression
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < VectorExpression<A,B,C>, DyadicOp, VectorExpression<D,E,F> > 
Dyadic(const VectorExpression<A,B,C> & a, 
       const VectorExpression<D,E,F> & b);

// Inverse of a matrix
inline MatrixExpression<> Inverse(const MatrixExpression<> &M);

// rotation matrix around vector V by an angle V.nrm()
inline MatrixExpression<> Rodrigues(const Vector &V);

// distance between two vectors
template<class A,int B,class C,class D,int E,class F>
inline DOUBLE dist(const VectorExpression<A,B,C> &a, 
		   const  VectorExpression<D,E,F> &b); 

// distance squared between two vectors
template<class A,int B,class C,class D,int E,class F>
inline DOUBLE dist2(const VectorExpression<A,B,C> &a, 
		    const  VectorExpression<D,E,F> &b);

// distance between two vectors with the first one shifted by s
template<class A,int B,class C,class D,int E,class F,class G,int H,class I>
inline DOUBLE distwithshift(const VectorExpression<A,B,C> & a, 
			    const VectorExpression<D,E,F> & b,
			    const VectorExpression<G,H,I> & s);

// output to stream

// handle difference between compilers cxx (compaq alpha) 
// and g++ (e.g. linux) regarding ostream definition
#ifdef __alpha
#include <iostream.h>
#define OSTREAM ostream
#else
#include <iostream>
#define OSTREAM std::ostream
#endif

// output for vector expressions
template<class A, int B, class C>
inline OSTREAM& operator<< ( OSTREAM& o, 
			     const VectorExpression<A,B,C> & t );

// output for matrix expressions
template<class A, int B, class C>
inline OSTREAM& operator<< ( OSTREAM& o, 
			     const MatrixExpression<A,B,C> & t );

//
// INLINE IMPLEMENTATION
//

// make sure a function for taking a square exists

#ifndef SQR
#ifdef _MYINLINEH_
#define SQR sqr
#else
inline DOUBLE vec3et_sqr(DOUBLE x) { return x*x; }
#define SQR vec3et_sqr
#endif
#endif

// make sure a SINCOS function is defined
#ifndef SINCOS
#define SINCOS1(t,s,c) sincos(t, s, c)
#define SINCOS2(t,s,c) sincos(s, c ,t)
#define SINCOS3(t,s,c) (*s=sin(t)),(*c=cos(t))
#ifdef __alpha
#define SINCOS(t,s,c) SINCOS3(t,s,c)
#else
#ifdef __GNUC__
#ifdef __DJGPP__
#define SINCOS(t,s,c) SINCOS2(t,s,c)
#else
#ifdef __USE_GNU
#define SINCOS(t,s,c) SINCOS1(t,s,c)
#else
#define SINCOS(t,s,c) SINCOS3(t,s,c)
#endif
#endif
#else
#define SINCOS(t,s,c) SINCOS3(t,s,c)
#endif
#endif
#endif

//
// member functions of the basic Vector type
//

// constructor
inline VectorExpression<>::VectorExpression (DOUBLE a, DOUBLE b, DOUBLE c)
  : x(a), y(b), z(c) {}

// constuctor from a vector expression
template<class A,int B,class C> 
inline VectorExpression<>::VectorExpression(const VectorExpression<A,B,C> &e) :
  x( e.eval<0>() ),
  y( e.eval<1>() ),
  z( e.eval<2>() )
{}

// assign zero to all elements
inline void VectorExpression<>::zero()
{
  x = y = z = 0;
}

// access to elements through parenthesis
// active: for assignment to ()
inline DOUBLE& VectorExpression<>:: operator()(int i) 
{ 
  return *(&x+i);
}
// passive
inline const DOUBLE& VectorExpression<>:: operator() (int i) const
{ 
  return *(&x+i);
}

// access to elements through parenthesis
// active: for assignment to ()
inline DOUBLE& VectorExpression<>:: operator[](int i) 
{ 
  return *(&x+i);
}
// passive
inline const DOUBLE& VectorExpression<>:: operator[](int i) const
{ 
  return *(&x+i);
}

// passive access to the elements through eval
template <> inline DOUBLE VectorExpression<>::eval<0>() const { return x; }
template <> inline DOUBLE VectorExpression<>::eval<1>() const { return y; }
template <> inline DOUBLE VectorExpression<>::eval<2>() const { return z; }

// norm squared
inline DOUBLE VectorExpression<>::nrm2() const 
{
   return x*x+y*y+z*z;
}

// norm 
inline DOUBLE VectorExpression<>::nrm() const
{
  DOUBLE biggest=x<0?-x:x;
  DOUBLE absval=y<0?-y:y;
  if (absval>biggest) biggest=absval;
  absval=z<0?-z:z;
  if (absval>biggest) biggest=absval;
  if (biggest!=0)
    return biggest*(DOUBLE)(sqrt(SQR(x/biggest)
				 +SQR(y/biggest)
				 +SQR(z/biggest)));
  else
    return 0;
}

// assignment from an expression
template<class A,int B,class C> 
inline VectorExpression<>& VectorExpression<>:: operator=(const VectorExpression<A,B,C> &e)
{
  x = e.eval<0>();
  y = e.eval<1>();
  z = e.eval<2>();
  return *this;
}

// addto-assignment from an expression
template<class A,int B,class C> 
inline VectorExpression<>& VectorExpression<>:: operator+=(const VectorExpression<A,B,C> &e)
{
  x += e.eval<0>();
  y += e.eval<1>();
  z += e.eval<2>();
  return *this;
}

// subtract-from assignment from an expression
template<class A,int B,class C>
inline VectorExpression<>& VectorExpression<>:: operator-=(const VectorExpression<A,B,C> &e)
{
  x -= e.eval<0>();
  y -= e.eval<1>();
  z -= e.eval<2>();
  return *this;
}

// multiply-by assignment from an expression
inline VectorExpression<>&  VectorExpression<>:: operator*=(const DOUBLE d)
{
  x *= d;
  y *= d;
  z *= d;
  return *this;
}

// divide-by assignment from an expression
inline VectorExpression<>&  VectorExpression<>:: operator/=(const DOUBLE d)
{
  x /= d;
  y /= d;
  z /= d;
  return *this;
}

//                                           
// member functions of the basic Matrix type 
//                                           

// constructor
inline MatrixExpression<>::MatrixExpression
  (DOUBLE a, DOUBLE b, DOUBLE c,
   DOUBLE d, DOUBLE e, DOUBLE f,
   DOUBLE g, DOUBLE h, DOUBLE i) :
    xx(a), xy(b), xz(c),
    yx(d), yy(e), yz(f),
    zx(g), zy(h), zz(i)
{}

// constructor from an MatrixExpression
template<class A,int B,class C>
  inline MatrixExpression<>::MatrixExpression(const MatrixExpression<A,B,C>&e):
  xx(e.eval<0,0>()), xy(e.eval<0,1>()), xz(e.eval<0,2>()),
  yx(e.eval<1,0>()), yy(e.eval<1,1>()), yz(e.eval<1,2>()),
  zx(e.eval<2,0>()), zy(e.eval<2,1>()), zz(e.eval<2,2>())
{}

// set all elements to zero
inline void MatrixExpression<>::zero()
{
  xx = xy = xz = yx = yy = yz = zx = zy = zz = 0;
}

// turn into an identity matrix
inline void MatrixExpression<>::one()
{
  xx = yy = zz = 1;
  xy = xz = yx = yz = zx = zy = 0;
}

// set the elements of a row equal to those of a vector
template<class A, int B, class C> inline 
void MatrixExpression<>::setRow(const int i, 
				const VectorExpression<A,B,C> &e)
{
  switch(i) {
  case 0: xx = e.eval<0>(); xy = e.eval<1>(); xz = e.eval<2>(); break;
  case 1: yx = e.eval<0>(); yy = e.eval<1>(); yz = e.eval<2>(); break;
  case 2: zx = e.eval<0>(); zy = e.eval<1>(); zz = e.eval<2>(); break;
  }
}

// set the elements of a column equal to those of a vector
template<class A, int B, class C> inline 
void MatrixExpression<>::setColumn(const int j,
				   const VectorExpression<A,B,C>&e)
{
  switch(j) {
  case 0: xx = e.eval<0>(); yx = e.eval<1>(); zx = e.eval<2>(); break;
  case 1: xy = e.eval<0>(); yy = e.eval<1>(); zy = e.eval<2>(); break;
  case 2: xz = e.eval<0>(); yz = e.eval<1>(); zz = e.eval<2>(); break;
  }
}

// access to elements through parenthesis
// active: for assignment to ()
inline DOUBLE& MatrixExpression<>:: operator()(int i, int j) 
{ 
  return *(&xx+3*i+j);
}
// passive:
inline const DOUBLE& MatrixExpression<>:: operator() (int i, int j) const
{ 
  return *(&xx+3*i+j);
}

// passive access to the elements through eval
template <> inline DOUBLE MatrixExpression<>::eval<0,0>() const { return xx; }
template <> inline DOUBLE MatrixExpression<>::eval<0,1>() const { return xy; }
template <> inline DOUBLE MatrixExpression<>::eval<0,2>() const { return xz; }
template <> inline DOUBLE MatrixExpression<>::eval<1,0>() const { return yx; }
template <> inline DOUBLE MatrixExpression<>::eval<1,1>() const { return yy; }
template <> inline DOUBLE MatrixExpression<>::eval<1,2>() const { return yz; }
template <> inline DOUBLE MatrixExpression<>::eval<2,0>() const { return zx; }
template <> inline DOUBLE MatrixExpression<>::eval<2,1>() const { return zy; }
template <> inline DOUBLE MatrixExpression<>::eval<2,2>() const { return zz; }

// trace
inline DOUBLE MatrixExpression<>::tr() const 
{ 
  return xx+yy+zz; 
}

// determininant
inline DOUBLE MatrixExpression<>::det() const
{
  return xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx);
}

// norm squared
inline DOUBLE MatrixExpression<>::nrm2() const 
{ 
  return SQR(xx)+SQR(xy)+SQR(xz)
        +SQR(yx)+SQR(yy)+SQR(yz)
        +SQR(zx)+SQR(zy)+SQR(zz);
}

// norm
inline DOUBLE MatrixExpression<>::nrm() const 
{
  DOUBLE absval;
  DOUBLE biggest=xx<0?-xx:xx;
  absval=xy<0?-xy:xy;
  if (absval>biggest) biggest=absval;
  absval=xz<0?-xz:xz;
  if (absval>biggest) biggest=absval;
  absval=yx<0?-yx:yx;
  if (absval>biggest) biggest=absval;
  absval=yy<0?-yy:yy;
  if (absval>biggest) biggest=absval;
  absval=yz<0?-yz:yz;
  if (absval>biggest) biggest=absval;
  absval=zx<0?-zx:zx;
  if (absval>biggest) biggest=absval;
  absval=zy<0?-zy:zy;
  if (absval>biggest) biggest=absval;
  absval=zz<0?-zz:zz;
  if (absval>biggest) biggest=absval;
  if (biggest!=0) 
    return biggest*(DOUBLE)(sqrt(SQR(xx/biggest)
			       +SQR(xy/biggest)
			       +SQR(xz/biggest)
			       +SQR(yx/biggest)
			       +SQR(yy/biggest)
			       +SQR(yz/biggest)
			       +SQR(zx/biggest)
			       +SQR(zy/biggest)
			       +SQR(zz/biggest)));
  else 
    return 0;
} 

// assignment from an expression
template<class A,int B,class C>
inline MatrixExpression<>& MatrixExpression<>:: operator=(const MatrixExpression<A,B,C> &e)
{
  xx = e.eval<0,0>(); xy = e.eval<0,1>(); xz = e.eval<0,2>();
  yx = e.eval<1,0>(); yy = e.eval<1,1>(); yz = e.eval<1,2>();
  zx = e.eval<2,0>(); zy = e.eval<2,1>(); zz = e.eval<2,2>();
  return *this;
}

// multiply-by assignment from an expression
inline MatrixExpression<>& MatrixExpression<>:: operator*=(const DOUBLE d)
{
  xx *= d;  xy *= d;  xz *= d;
  yx *= d;  yy *= d;  yz *= d;
  zx *= d;  zy *= d;  zz *= d;
  return *this;
}
// divide-by assignment from an expression
inline MatrixExpression<>& MatrixExpression<>:: operator/=(const DOUBLE d)
{
  xx /= d;  xy /= d;  xz /= d;
  yx /= d;  yy /= d;  yz /= d;
  zx /= d;  zy /= d;  zz /= d;
  return *this;
}

// add-to assignment from an expression
template<class A,int B,class C>
inline MatrixExpression<>& MatrixExpression<>:: operator+=(const MatrixExpression<A,B,C> &e)
{
  xx += e.eval<0,0>(); xy += e.eval<0,1>(); xz += e.eval<0,2>();
  yx += e.eval<1,0>(); yy += e.eval<1,1>(); yz += e.eval<1,2>();
  zx += e.eval<2,0>(); zy += e.eval<2,1>(); zz += e.eval<2,2>();
  return *this;
}
  
// subtract-from assignment from an expression
template<class A,int B,class C>
inline MatrixExpression<>& MatrixExpression<>:: operator-=(const MatrixExpression<A,B,C> &e)
{
  xx -= e.eval<0,0>(); xy -= e.eval<0,1>(); xz -= e.eval<0,2>();
  yx -= e.eval<1,0>(); yy -= e.eval<1,1>(); yz -= e.eval<1,2>();
  zx -= e.eval<2,0>(); zy -= e.eval<2,1>(); zz -= e.eval<2,2>();
  return *this;
}

// multiply-by-matrix assignment from an expression
template<class A,int B,class C>
inline MatrixExpression<>& MatrixExpression<>::operator*=(const MatrixExpression<A,B,C> & m)
{
  MatrixExpression<> rhs(m);
  DOUBLE x = xx, y = xy;
  xx *= rhs.xx; xx += y*rhs.yx+xz*rhs.zx;
  xy *= rhs.yy; xy += x*rhs.xy+xz*rhs.zy;
  xz *= rhs.zz; xz += x*rhs.xz+ y*rhs.yz;
  x = yx; y = yy;
  yx *= rhs.xx; yx += y*rhs.yx+yz*rhs.zx;
  yy *= rhs.yy; yy += x*rhs.xy+yz*rhs.zy;
  yz *= rhs.zz; yz += x*rhs.xz+ y*rhs.yz;
  x = zx; y = zy;
  zx *= rhs.xx; zx += y*rhs.yx+zz*rhs.zx;
  zy *= rhs.yy; zy += x*rhs.xy+zz*rhs.zy;
  zz *= rhs.zz; zz += x*rhs.xz+ y*rhs.yz;
  return *this;
}


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

// to access the elements of a matrix expression
#define MATPARENTHESIS operator()(const int _i, const int _j) const \
{ \
  switch(_i){ \
    case 0: switch(_j){ \
              case 0: return eval<0,0>(); \
              case 1: return eval<0,1>(); \
              case 2: return eval<0,2>(); \
              default: return 0; \
            } \
    case 1: switch(_j){ \
              case 0: return eval<1,0>(); \
              case 1: return eval<1,1>(); \
              case 2: return eval<1,2>(); \
              default: return 0; \
	    } \
    case 2: switch(_j){ \
              case 0: return eval<2,0>(); \
              case 1: return eval<2,1>(); \
              case 2: return eval<2,2>(); \
              default: return 0; \
            } \
    default: return 0; \
  } \
}

// to get the norm of a vector
#define VECNRM nrm() const\
{\
  DOUBLE x=eval<0>(); \
  DOUBLE y=eval<1>(); \
  DOUBLE z=eval<2>(); \
  DOUBLE biggest=x<0?-x:x; \
  DOUBLE absval=y<0?-y:y; \
  if (absval>biggest) biggest=absval; \
  absval=z<0?-z:z; \
  if (absval>biggest) biggest=absval; \
  if (biggest!=0) \
    return biggest*(DOUBLE)(sqrt(SQR(x/biggest) \
			       +SQR(y/biggest) \
			       +SQR(z/biggest))); \
  else \
    return 0; \
}

// to get the norm of a matrix
#define MATNRM nrm() const \
{ \
  DOUBLE xx=eval<0,0>(), xy=eval<0,1>(), xz=eval<0,2>(); \
  DOUBLE yx=eval<1,0>(), yy=eval<1,1>(), yz=eval<1,2>(); \
  DOUBLE zx=eval<2,0>(), zy=eval<2,1>(), zz=eval<2,2>(); \
  DOUBLE absval;\
  DOUBLE biggest=xx<0?-xx:xx;\
  absval=xy<0?-xy:xy;\
  if (absval>biggest) biggest=absval;\
  absval=xz<0?-xz:xz;\
  if (absval>biggest) biggest=absval;\
  absval=yx<0?-yx:yx;\
  if (absval>biggest) biggest=absval;\
  absval=yy<0?-yy:yy;\
  if (absval>biggest) biggest=absval;\
  absval=yz<0?-yz:yz;\
  if (absval>biggest) biggest=absval;\
  absval=zx<0?-zx:zx;\
  if (absval>biggest) biggest=absval;\
  absval=zy<0?-zy:zy;\
  if (absval>biggest) biggest=absval;\
  absval=zz<0?-zz:zz;\
  if (absval>biggest) biggest=absval;\
  if (biggest!=0) \
    return biggest*(DOUBLE)(sqrt(SQR(xx/biggest)\
			       +SQR(xy/biggest)\
			       +SQR(xz/biggest)\
			       +SQR(yx/biggest)\
			       +SQR(yy/biggest)\
			       +SQR(yz/biggest)\
			       +SQR(zx/biggest)\
			       +SQR(zy/biggest)\
			       +SQR(zz/biggest))); \
  else                                    \
    return 0;                           \
} 


// to get the norm squared of a vector expression
#define VECNRM2 nrm2() const \
{\
   return SQR(eval<0>()) + SQR(eval<1>()) + SQR(eval<2>()); \
}

// to get the norm squared of a matrix expression
#define MATNRM2 nrm2() const \
{ \
  return SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>())  \
        +SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>())  \
        +SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>()); \
}

// to get the trace of a matrix expression
#define MATTR tr() const \
{ \
  return eval<0,0>() + eval<1,1>() + eval<2,2>(); \
}

// to get the determinant of a matrix expression
#define MATDET det() const \
{ \
  DOUBLE xx=eval<0,0>(), xy=eval<0,1>(), xz=eval<0,2>(); \
  DOUBLE yx=eval<1,0>(), yy=eval<1,1>(), yz=eval<1,2>(); \
  DOUBLE zx=eval<2,0>(), zy=eval<2,1>(), zz=eval<2,2>(); \
  return xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx); \
}

// to get the i-th row of a matrix expression
// note: the class of the matrix expression has to be defined in CLASS
#define MATROW VectorExpression<CLASS,RowOp,Base> row(int i) const \
{ \
  return VectorExpression<CLASS,RowOp,Base>( *this, i ); \
}

// to get the j-th column of a matrix expression
// note: the class of the matrix expression has to be defined in CLASS
#define MATCOLUMN VectorExpression<CLASS,ColumnOp,Base> column(int j) const \
{ \
  return VectorExpression<CLASS,ColumnOp,Base>( *this, j ); \
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
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to
  template<int I> inline DOUBLE eval() const
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
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to
  template <int I> inline DOUBLE eval() const
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
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES
 
  // define what this operation evaluates to
  template <const int I> inline DOUBLE eval() const
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


// Expression template class for '*' between a Vector and a DOUBLE
//
// A and C are classes, B is an Operator
// if B=NoOp and A=C=Base, the VectorExpression is actually a Vector
template<class A,int B,class C>
class VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
{
 public:

  // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
  inline DOUBLE VECNRM2  
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES
 
  // define what this operation evaluates to
  template <int I> inline DOUBLE eval() const
  { 
    return a->eval<I>() * b; 
  }

  // optimize a second multiplication with a DOUBLE
  inline VectorExpression & operator* (DOUBLE c)
  {
    b *= c;
    return *this;
  }

  // optimize a subsequent division by a DOUBLE
  inline VectorExpression & operator/ (DOUBLE c)
  {
    b /= c;
    return *this;
  }

  // constructor
  inline VectorExpression(const VectorExpression<A,B,C> & _a, 
			  DOUBLE _b) 
    : a(&_a), b(_b) {}

 private:
  // pointer to the subexpression 
  const VectorExpression<A,B,C>* a;
  // the DOUBLE to multiply with
  DOUBLE b;

  // be-friend multiplication operators that optimize further DOUBLE
  // multiplication
  template<class D,int E,class F> 
  friend VectorExpression < VectorExpression<D,E,F>, TimesOp, DOUBLE >&
  operator* (DOUBLE b, VectorExpression < VectorExpression<D,E,F>, TimesOp, DOUBLE > &a);

  template<class D,int E,class F> 
  friend VectorExpression < VectorExpression<D,E,F>, TimesOp, DOUBLE >&
  operator/ (DOUBLE b, VectorExpression < VectorExpression<D,E,F>, TimesOp, DOUBLE > &a);

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
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to
  template <int I> inline DOUBLE eval() const
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


// Expression template class for '+' between two MatrixExpressions
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Matrix
// if E=NoOp and D=F=Base, the second sub-expression is actually a Matrix
template<class A,int B,class C,class D,int E,class F>
class MatrixExpression < MatrixExpression<A,B,C>, PlusOp, MatrixExpression<D,E,F> >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(), element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS 

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,PlusOp,MatrixExpression<D,E,F> >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  { 
    return a->eval<I,J>() + b->eval<I,J>(); 
  }

  // constructor
  inline 
  MatrixExpression(const MatrixExpression<A,B,C> & _a, 
		   const MatrixExpression<D,E,F> & _b)
    : a(&_a), b(&_b) {}

 private:
  // store pointers to the sub-expressions
  const MatrixExpression<A,B,C> * a;
  const MatrixExpression<D,E,F> * b;
};


// Expression template class for '-' between two MatrixExpressions
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Matrix
// if E=NoOp and D=F=Base, the second sub-expression is actually a Matrix
template<class A,int B,class C,class D,int E,class F>
class MatrixExpression < MatrixExpression<A,B,C>, MinusOp, MatrixExpression<D,E,F> >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(), and
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,MinusOp,MatrixExpression<D,E,F> >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  { 
    return a->eval<I,J>() - b->eval<I,J>(); 
  }

  // constructor
  inline 
  MatrixExpression(const MatrixExpression<A,B,C>& _a, 
		   const MatrixExpression<D,E,F>& _b)
    : a(&_a), b(&_b) {}

 private:
  // store pointers to the sub-expressions
  const MatrixExpression<A,B,C> * a;
  const MatrixExpression<D,E,F> * b;
};


// Expression template class for '*' between a MatrixExpression and a DOUBLE
//
// A and C are classes, B is an operator
// if B=NoOp and A=C=Base, the sub-expression is actually a Matrix
template<class A,int B,class C>
class MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
{
 public:
  
  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,TimesOp,DOUBLE >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // constructor
  inline 
  MatrixExpression(const MatrixExpression<A,B,C> & _a, 
		   DOUBLE _b) 
    : a(&_a), b(_b) {}

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  { return a->eval<I,J>() * b; }

  inline MatrixExpression & operator* (DOUBLE c)
  {
    b *= c;
    return *this;
  }

  inline MatrixExpression & operator/ (DOUBLE c)
  {
    b /= c;
    return *this;
  }

 private:
  // store the DOUBLE and a pointer to the sub-MatrixExpression
  const MatrixExpression<A,B,C>* a;
  DOUBLE b;

  template<class D,int E,class F> 
  friend MatrixExpression < MatrixExpression<D,E,F>, TimesOp, DOUBLE >&
  operator* (DOUBLE b, MatrixExpression < MatrixExpression<D,E,F>, TimesOp, DOUBLE > &a);

  template<class D,int E,class F> 
  friend MatrixExpression < MatrixExpression<D,E,F>, TimesOp, DOUBLE >&
  operator/ (DOUBLE b, MatrixExpression < MatrixExpression<D,E,F>, TimesOp, DOUBLE > &a);
};


// Expression template class for '*' between two MatrixExpressions
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Matrix
// if E=NoOp and D=F=Base, the second sub-expression is actually a Matrix
template<class A,int B,class C,class D,int E,class F>
class MatrixExpression < MatrixExpression<A,B,C>, TimesOp, MatrixExpression<D,E,F> >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,TimesOp,MatrixExpression<D,E,F> >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  {
    return  a->eval<I,0>() * b->eval<0,J>()
          + a->eval<I,1>() * b->eval<1,J>()
          + a->eval<I,2>() * b->eval<2,J>();
  }

  // constructor
  inline MatrixExpression(const MatrixExpression<A,B,C>&_a, 
			  const MatrixExpression<D,E,F>&_b) 
    : a(&_a), b(&_b) {}

 private:
  // store pointers to the sub-expressions
  const MatrixExpression<A,B,C>* a;
  const MatrixExpression<D,E,F>* b;
};

// Expression template class for unary '-' acting on a MatrixExpression 
//
// A and C are classes, B is an operator
// if B=NoOp and A=C=Base, the sub-expression is actually a Matrix
template<class A,int B,class C>
class MatrixExpression < MatrixExpression<A,B,C>, NegativeOp, Base >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< MatrixExpression<A,B,C>, NegativeOp, Base >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  {
    return  - a->eval<I,J>;
  }

  // constructor
  inline MatrixExpression(const MatrixExpression<A,B,C> &_a) 
    : a(&_a) {}

 private:
  // store pointer to the sub-expression
  const MatrixExpression<A,B,C> * a;
};


// Expression template class for transpose function 
//
// A and C are classes, B is an operator
// if B=NoOp and A=C=Base, the sub-expression is actually a Matrix
template<class A,int B,class C>
class MatrixExpression < MatrixExpression<A,B,C>, TransposeOp, Base >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< MatrixExpression<A,B,C>, TransposeOp, Base >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  {
    return  a->eval<J,I>();
  }

  // constructor
  inline 
  MatrixExpression(const MatrixExpression<A,B,C> &_a) 
    : a(&_a) {}

 private:
  // store pointer to the sub-expression
  const MatrixExpression<A,B,C> * a;
};


// Expression template class for Dyadic operation between two VectorExpressions
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Vector
// if E=NoOp and D=F=Base, the second sub-expression is actually a Vector
template<class A,int B,class C,class D,int E,class F>
class MatrixExpression < VectorExpression<A,B,C>, DyadicOp, VectorExpression<D,E,F> >
{
 public:

  // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
  // element access (*this)(i,j)
  inline DOUBLE MATNRM2
  inline DOUBLE MATNRM
  inline DOUBLE MATDET
  inline DOUBLE MATTR
  inline DOUBLE MATPARENTHESIS

  // define row and column access as if they were Vectors
  // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< VectorExpression<A,B,C>, DyadicOp, VectorExpression<D,E,F> >
  inline MATROW
  inline MATCOLUMN
#undef CLASS

  // define what this operation evaluates to  
  template <int I, int J> inline DOUBLE eval() const
  {
    return a->eval<I>() * b->eval<J>();
  }

  // constructor
  inline 
  MatrixExpression(const VectorExpression<A,B,C> & _a, 
		   const VectorExpression<D,E,F> & _b) 
    : a(&_a), b(&_b) {}

 private:
  // store pointers to the sub-expressions
  const VectorExpression<A,B,C> * a;
  const VectorExpression<D,E,F> * b;
};


// Expression template class for row operation acting on a MatrixExpression 
//
// A and C are classes, B is an operator
// if B=NoOp and A=C=Base, the sub-expression is actually a Matrix
template<class A,int B,class C>
class VectorExpression < MatrixExpression<A,B,C>, RowOp, Base >
{
 public:

  // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to  
  template <const int J> inline DOUBLE eval() const
  {
    switch(i) {
      case 0: return a->eval<0,J>();
      case 1: return a->eval<1,J>();
      case 2: return a->eval<2,J>();
      default: 
	return 0; // should never happen
    }
  }

  // constructor
  inline 
  VectorExpression(const MatrixExpression<A,B,C> & _a,
		   int _i) 
    : a(&_a), i(_i) {}

 private:
  // store pointer to the sub-expression and the index i
  const MatrixExpression<A,B,C> * a;
  const int i;
};


// Expression template class for column operation acting on a MatrixExpression 
//
// A and C are classes, B is an operator
// if B=NoOp and A=C=Base, the sub-expression is actually a Matrix
template<class A,int B,class C>
class VectorExpression < MatrixExpression<A,B,C>, ColumnOp, Base>
{
 public:

  // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to  
  template <const int I> inline DOUBLE eval() const
  {
    switch(j) {
      case 0: return a->eval<I,0>();
      case 1: return a->eval<I,1>();
      case 2: return a->eval<I,2>();
      default: 
	return 0; // should never happen
    }
  }

  // constructor
  inline 
  VectorExpression(const MatrixExpression<A,B,C> & _a, 
		   int _j) 
    :  a(&_a), j(_j) {}

 private:
  // store pointer to the sub-expression and the index j
  const MatrixExpression<A,B,C> * a;
  const int j;
};


// Expression template class for '*' between a MatrixExpression and a
// VectorExpression
//
// A, C, D, and F are classes, B and E are Operators
// if B=NoOp and A=C=Base, the first sub-expression is actually a Matrix
// if E=NoOp and D=F=Base, the second sub-expression is actually a Vector
template<class A,int B,class C,class D,int E,class F>
class VectorExpression < MatrixExpression<A,B,C>, TimesOp, VectorExpression <D,E,F> >
{
 public:

  // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
  inline DOUBLE VECNRM2
  inline DOUBLE VECNRM
  inline DOUBLE VECPARENTHESES

  // define what this operation evaluates to  
  template <int I> inline DOUBLE eval() const
  {
    return  a->eval<I,0>() * b->eval<0>()
          + a->eval<I,1>() * b->eval<1>()
          + a->eval<I,2>() * b->eval<2>();
  }

  // constructor
  inline 
  VectorExpression(const MatrixExpression<A,B,C> &_a, 
		   const VectorExpression <D,E,F> &_b) 
    : a(&_a), b(&_b) {}

 private:
  // store pointers to the sub-expressions
  const MatrixExpression<A,B,C> * a;
  const VectorExpression<D,E,F> * b;
};

// helper class to implemement a comma list assignment
class CommaOp 
{
 public:

  // the comma operatorset one element and moves the pointer to the next one
  inline CommaOp& operator,(DOUBLE c) 
  { 
    // make sure the vector or matrix is not completely filled yet
    if (ptr <= end) 
      *ptr++ = c; 
    return *this; 
  }

  // destructor: fill whatever elements remain with zeros
  inline ~CommaOp() 
  { 
    while (ptr <= end) 
      *ptr++ = 0; 
  }

 private:
  // store pointers to the next element to be filled and the last
  // element fillable.
  DOUBLE * ptr;
  DOUBLE * const end; 

  // constructor: set pointer to the element to be filled next and to
  // the last element
  inline CommaOp(DOUBLE & a,DOUBLE & b) 
    : ptr(&a), end(&b) {}

  // this constructor is private, and so that the CommaOp class can
  // only be used by the friend member functions operator= of
  // VectorExpression<> and MatrixExpression<>
  friend CommaOp VectorExpression<>::operator=(DOUBLE c);
  friend CommaOp MatrixExpression<>::operator=(DOUBLE c);
};

// assignment operatorof VectorExpression<> that triggers a comma
// separated initializer list
inline CommaOp VectorExpression<>::operator=(DOUBLE c)
{
  x = c;
  return CommaOp(y, z);
}

// assignment operatorof MatrixExpression<> that triggers a comma
// separated initializer list
inline CommaOp MatrixExpression<>::operator=(DOUBLE c)
{
  xx = c;
  return CommaOp(xy, zz);
}


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

// DOUBLE * Vector
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator* (DOUBLE b, 
	   const VectorExpression<A,B,C> & a)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>(a, b);
}

// Vector * DOUBLE
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator* (const VectorExpression<A,B,C> & a, 
	   DOUBLE b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>(a, b);
}

// Vector / DOUBLE
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>
operator/ (const VectorExpression<A,B,C> & a, 
	   DOUBLE b)
{ 
  return VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>(a, 1.0/b);
}

// Vector * DOUBLE * DOUBLE
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>&
operator*(DOUBLE b,VectorExpression<VectorExpression<A,B,C>,TimesOp,DOUBLE>&a)
{ 
  a.b *= b; 
  return a; 
}

// Vector * DOUBLE / DOUBLE
template<class A,int B,class C> 
inline VectorExpression<VectorExpression<A,B,C>, TimesOp, DOUBLE>&
operator/(DOUBLE b,VectorExpression<VectorExpression<A,B,C>,TimesOp,DOUBLE>&a)
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
inline DOUBLE 
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
inline DOUBLE 
operator*(const VectorExpression<A,B,C> & a, const VectorExpression<D,E,F> & b)
{
  return 
    a.template eval<0>()*b.template eval<0>() 
    + a.template eval<1>()*b.template eval<1>() 
    + a.template eval<2>()*b.template eval<2>();
}

// Matrix + Matrix
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, PlusOp, MatrixExpression<D,E,F> >
operator+ (const MatrixExpression<A,B,C> & a, 
	   const MatrixExpression<D,E,F> & b)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, PlusOp, MatrixExpression<D,E,F> >(a,b);
}

// Matrix - Matrix
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, MinusOp, MatrixExpression<D,E,F> >
operator- (const MatrixExpression<A,B,C> & a, 
	   const MatrixExpression<D,E,F> & b)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, MinusOp, MatrixExpression<D,E,F> >(a, b);
}

// DOUBLE * Matrix
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator* (DOUBLE b, 
	   const MatrixExpression<A,B,C> & a)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >(a, b);
}

// Matrix * DOUBLE
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator* (const MatrixExpression<A,B,C> & a, 
	   DOUBLE b)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >(a, b);
}

// Matrix / DOUBLE
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >
operator/ (const MatrixExpression<A,B,C> & a, 
	   DOUBLE b)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >(a, 1.0/b);
}

// Matrix * DOUBLE * DOUBLE
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >&
operator* (DOUBLE b, MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE > &a)
{ 
  a.b *= b; 
  return a; 
}

// Matrix * DOUBLE / DOUBLE
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE >&
operator/ (DOUBLE b, MatrixExpression < MatrixExpression<A,B,C>, TimesOp, DOUBLE > &a)
{ 
  a.b /= b; 
  return a; 
}

// Matrix * Matrix
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < MatrixExpression<A,B,C>, TimesOp, MatrixExpression<D,E,F> >
operator*(const MatrixExpression<A,B,C> & a, 
	  const MatrixExpression<D,E,F> & b)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, TimesOp, MatrixExpression<D,E,F> > (a, b);
}

// -Matrix 
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, NegativeOp, Base >
operator- (const MatrixExpression<A,B,C> & a)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, NegativeOp, Base > (a);
}

// Matrix * Vector
template<class A,int B,class C,class D,int E,class F>
inline VectorExpression < MatrixExpression<A,B,C>, TimesOp, VectorExpression<D,E,F> >
operator*(const MatrixExpression<A,B,C> & a, 
	  const VectorExpression<D,E,F> & b)
{ 
  return VectorExpression < MatrixExpression<A,B,C>, TimesOp, VectorExpression<D,E,F> > (a, b);
}

//
// Implementation of the row and column member functions of MatrixExpression
//

inline VectorExpression<MatrixExpression<>,RowOp,Base> 
MatrixExpression<>::row(int i)
{
  return VectorExpression<MatrixExpression<>,RowOp,Base>(*this,i);
}

inline VectorExpression<MatrixExpression<>,ColumnOp,Base> 
MatrixExpression<>::column(int j)
{
  return VectorExpression<MatrixExpression<>,ColumnOp,Base>(*this,j);
}

// Can now implement Gramm-Schmidt orthogonalization of the rows of the matrix
inline void Matrix::reorthogonalize()
{
  DOUBLE z;
  int num=10;
  while ( fabs(det()-1) > 1E-16 && --num) {
    z = 1/row(0).nrm();   
    xx *= z;    
    xy *= z;    
    xz *= z;    
    z = xx*yx + xy*yy + xz*yz;          
    yx -= z*xx; 
    yy -= z*xy; 
    yz -= z*xz; 
    z = 1/row(1).nrm();                
    yx *= z;  
    yy *= z;  
    yz *= z;        
    zx = xy*yz - xz*yy;                 
    zy = xz*yx - xx*yz;
    zz = xx*yy - xy*yx;
    z = 1/row(2).nrm();                
    zx *= z;  
    zy *= z;  
    zz *= z;        
  }
}

//
// Non-member functions
//

// return |a-b|
template<class A,int B,class C,class D,int E,class F>
inline DOUBLE dist(const VectorExpression<A,B,C> & a, 
		   const VectorExpression<D,E,F> & b)
{
  DOUBLE x = fabs(a.template eval<0>()-b.template eval<0>());
  DOUBLE y = fabs(a.template eval<1>()-b.template eval<1>());
  DOUBLE z = fabs(a.template eval<2>()-b.template eval<2>());
  DOUBLE biggest=x;
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
inline DOUBLE dist2(const VectorExpression<A,B,C> & a, 
		    const VectorExpression<D,E,F> & b)
{
  DOUBLE d = a.template eval<0>() - b.template eval<0>();
  DOUBLE c = d * d;
  d = a.template eval<1>() - b.template eval<1>();
  c += d * d;
  d = a.template eval<2>() - b.template eval<2>();
  return c + d * d;
}

// return |a+s-b|
template<class A,int B,class C,class D,int E,class F,class G,int H,class I>
inline DOUBLE distwithshift(const VectorExpression<A,B,C> &a, 
			    const VectorExpression<D,E,F> & b, 
			    const VectorExpression<G,H,I>& s)
{
 DOUBLE x=fabs(s.template eval<0>()+a.template eval<0>()-b.template eval<0>());
 DOUBLE y=fabs(s.template eval<1>()+a.template eval<1>()-b.template eval<1>());
 DOUBLE z=fabs(s.template eval<2>()+a.template eval<2>()-b.template eval<2>());
 DOUBLE biggest=x;
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

// the inverse of a matrix
inline MatrixExpression<> Inverse(const MatrixExpression<> &M) 
{
  DOUBLE d = 1/M.det();
  return MatrixExpression<>(
    (M.yy*M.zz-M.yz*M.zy)*d, -(M.xy*M.zz-M.xz*M.zy)*d,  (M.xy*M.yz-M.xz*M.yy)*d,
   -(M.yx*M.zz-M.yz*M.zx)*d,  (M.xx*M.zz-M.xz*M.zx)*d, -(M.xx*M.yz-M.xz*M.yx)*d,
    (M.yx*M.zy-M.yy*M.zx)*d, -(M.xx*M.zy-M.xy*M.zx)*d,  (M.xx*M.yy-M.xy*M.yx)*d);
}

// rotation matrix built using the Rodrigues formula
inline MatrixExpression<> Rodrigues(const Vector &V) 
{
  DOUBLE theta = V.nrm();
  if (theta != 0) {
    DOUBLE s, c;
    SINCOS(theta,&s,&c);
    DOUBLE inrm = 1/theta;
    DOUBLE wx = V.x*inrm;
    DOUBLE wy = V.y*inrm;
    DOUBLE wz = V.z*inrm;
    DOUBLE oneminusc = 1-c;
    DOUBLE wxwy1mc = wx*wy*oneminusc;
    DOUBLE wxwz1mc = wx*wz*oneminusc;
    DOUBLE wywz1mc = wy*wz*oneminusc;
    DOUBLE wxs = wx*s;
    DOUBLE wys = wy*s;
    DOUBLE wzs = wz*s;
    return MatrixExpression<>(c+wx*wx*oneminusc,wxwy1mc-wzs,      wxwz1mc+wys,
			      wxwy1mc+wzs,      c+wy*wy*oneminusc,wywz1mc-wxs,
			      wxwz1mc-wys,     wywz1mc+wxs, c+wz*wz*oneminusc);
  } else {
    return MatrixExpression<>(1,0,0,0,1,0,0,0,1);
  }
}

// transpose of a matrix
template<class A,int B,class C>
inline MatrixExpression < MatrixExpression<A,B,C>, TransposeOp, Base >
Transpose (const MatrixExpression<A,B,C> & a)
{ 
  return MatrixExpression < MatrixExpression<A,B,C>, TransposeOp, Base > (a);
}

// dyadic product of two vectors
template<class A,int B,class C,class D,int E,class F>
inline MatrixExpression < VectorExpression<A,B,C>, DyadicOp, VectorExpression<D,E,F> >
Dyadic (const VectorExpression<A,B,C> & a, 
	const VectorExpression<D,E,F> & b)
{ 
  return MatrixExpression < VectorExpression<A,B,C>, DyadicOp, VectorExpression<D,E,F> > (a,b);
}

//
// Output
//

// vectors
template<class A, int B, class C>
inline OSTREAM& operator<< ( OSTREAM & o, 
			     const VectorExpression<A,B,C> & t )
{
  return o << t.template eval<0>() << ' ' 
	   << t.template eval<1>() << ' ' 
	   << t.template eval<2>();
}

// matrices
template<class A,int B,class C>
OSTREAM& operator<< ( OSTREAM & o, 
		      const MatrixExpression<A,B,C> & t )
{
  return o << '\n'<< t.template eval<0,0>() 
	   << ' ' << t.template eval<0,1>() 
	   << ' ' << t.template eval<0,2>()
	   << '\n'<< t.template eval<1,0>() 
	   << ' ' << t.template eval<1,1>() 
	   << ' ' << t.template eval<1,2>()
	   << '\n'<< t.template eval<2,0>() 
	   << ' ' << t.template eval<2,1>() 
	   << ' ' << t.template eval<2,2>() 
	   << '\n';
}

#undef OSTREAM
#undef VECNRM
#undef VECNRM2
#undef VECPARENTHESIS
#undef MATNRM
#undef MATNRM2
#undef MATPARENTHESIS
#undef MATTR
#undef MATDET
#undef MATROW
#undef MATCOLUMN


#endif

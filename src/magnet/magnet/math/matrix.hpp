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
#include <magnet/math/vector.hpp>
#include <magnet/math/detail/eigenval.hpp>

namespace magnet {
  namespace math {
    // define default MatrixExpression, which is the Matrix type
    template <>
    class MatrixExpression<>
    {
    public:
      double xx, xy, xz,                         // elements of the matrix
	yx, yy, yz,
	zx, zy, zz;

      inline MatrixExpression<>() {}             // default constructor

      inline MatrixExpression
      (double a, double b, double c,
       double d, double e, double f,
       double g, double h, double i) :
	xx(a), xy(b), xz(c),
	yx(d), yy(e), yz(f),
	zx(g), zy(h), zz(i)
      {}
      
      std::pair<std::tr1::array<Vector, 3>, std::tr1::array<double, 3> > 
      symmetric_eigen_decomposition() const
      {
#ifdef MAGNET_DEBUG
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    if (operator()(i,j) != operator()(i,j))
	      M_throw() << "Cannot perform an eigen decomposition of a matrix which is not symmetric using this function!";
#endif

	double V[3][3];
	double d[3];
	double e[3];
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    V[i][j] = operator()(i,j);
	
	detail::tred2(V, d, e);
	detail::tql2(V, d, e);

	std::tr1::array<double, 3> eigenvals = {{d[0], d[1], d[2]}};
	std::tr1::array<Vector, 3> eigenvecs 
	  = {{Vector(V[0][0], V[1][0], V[2][0]),
	      Vector(V[0][1], V[1][1], V[2][1]),
	      Vector(V[0][2], V[1][2], V[2][2])}};
	//Now output the eigenvectors and eigenvalues
	return std::make_pair(eigenvecs, eigenvals);
      }


      template<class A,int B,class C>
      inline MatrixExpression(const MatrixExpression<A,B,C>&e):
	xx(e.eval<0,0>()), xy(e.eval<0,1>()), xz(e.eval<0,2>()),
	yx(e.eval<1,0>()), yy(e.eval<1,1>()), yz(e.eval<1,2>()),
	zx(e.eval<2,0>()), zy(e.eval<2,1>()), zz(e.eval<2,2>())
      {}

      // set all elements to zero
      inline void zero()
      { xx = xy = xz = yx = yy = yz = zx = zy = zz = 0; }

      inline void one() // turn into an identity matrix
      { xx = yy = zz = 1; xy = xz = yx = yz = zx = zy = 0; }

      inline void reorthogonalize();             // Gramm-Schmidt
    
      inline double nrm2() const // norm squared
      { 
	return SQR(xx)+SQR(xy)+SQR(xz)
	  +SQR(yx)+SQR(yy)+SQR(yz)
	  +SQR(zx)+SQR(zy)+SQR(zz);
      }

      inline double tr() const { return xx+yy+zz; }  // trace

      // determininant
      inline double det() const
      { return xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx); }
  

      // norm
      inline double nrm() const 
      {
	double biggestRow1=std::max(std::max(std::abs(xx), std::abs(xy)), std::abs(xz));
	double biggestRow2=std::max(std::max(std::abs(yx), std::abs(yy)), std::abs(yz));
	double biggestRow3=std::max(std::max(std::abs(zx), std::abs(zy)), std::abs(zz));
	double biggest = std::max(std::max(biggestRow1, biggestRow2), biggestRow3);
	if (biggest!=0) 
	  return biggest*(double)(std::sqrt(SQR(xx/biggest)
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

      static inline Matrix identity() { return Matrix(1,0,0,0,1,0,0,0,1); }

      // access to elements through parenthesis
      inline double& operator()(int i, int j) { return *(&xx+3*i+j); }
      inline const double& operator() (int i, int j) const { return *(&xx+3*i+j); }  

      inline double& operator() (int i) { return *(&xx+i); }  
      inline const double& operator() (int i) const { return *(&xx+i); }  

      inline bool operator==(const MatrixExpression& omat) const
      {
	for (size_t i(0); i < 9; ++i)
	  if (operator()(i) != omat(i))
	    return false;

	return true;
      }

      inline VectorExpression<MatrixExpression<>,ops::RowOp,Base> row(int i);// return i-th row   
  
      inline VectorExpression<MatrixExpression<>,ops::ColumnOp,Base> column(int j); // return j-th column

      // set the elements of a row equal to those of a vector
      template<class A, int B, class C> inline 
      void setRow(const int i, const VectorExpression<A,B,C> &e)
      {
	switch(i) {
	case 0: xx = e.eval<0>(); xy = e.eval<1>(); xz = e.eval<2>(); break;
	case 1: yx = e.eval<0>(); yy = e.eval<1>(); yz = e.eval<2>(); break;
	case 2: zx = e.eval<0>(); zy = e.eval<1>(); zz = e.eval<2>(); break;
	}
      }

      // set the elements of a column equal to those of a vector
      template<class A, int B, class C> inline 
      void setColumn(const int j, const VectorExpression<A,B,C>&e)
      {
	switch(j) {
	case 0: xx = e.eval<0>(); yx = e.eval<1>(); zx = e.eval<2>(); break;
	case 1: xy = e.eval<0>(); yy = e.eval<1>(); zy = e.eval<2>(); break;
	case 2: xz = e.eval<0>(); yz = e.eval<1>(); zz = e.eval<2>(); break;
	}
      }

      inline MatrixExpression<>& operator*=(const double d);  // multiply by number
      inline MatrixExpression<>& operator/=(const double d);  // divide by number

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

      template<int I,int J> inline double eval() const; //evaluate components
    };

    // passive access to the elements through eval
    template <> inline double MatrixExpression<>::eval<0,0>() const { return xx; }
    template <> inline double MatrixExpression<>::eval<0,1>() const { return xy; }
    template <> inline double MatrixExpression<>::eval<0,2>() const { return xz; }
    template <> inline double MatrixExpression<>::eval<1,0>() const { return yx; }
    template <> inline double MatrixExpression<>::eval<1,1>() const { return yy; }
    template <> inline double MatrixExpression<>::eval<1,2>() const { return yz; }
    template <> inline double MatrixExpression<>::eval<2,0>() const { return zx; }
    template <> inline double MatrixExpression<>::eval<2,1>() const { return zy; }
    template <> inline double MatrixExpression<>::eval<2,2>() const { return zz; }


    // operator+ takes two matrix expressions and returns an appropriate
    // matrix expression
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::PlusOp, MatrixExpression<D,E,F> >
    operator+(const MatrixExpression<A,B,C> & a, 
	      const MatrixExpression<D,E,F> & b);

    // operator- takes two matrix expressions and returns an appropriate
    // matrix expression
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::MinusOp, MatrixExpression<D,E,F> >
    operator- (const MatrixExpression<A,B,C> & a,
	       const MatrixExpression<D,E,F> & b);

    // operator* takes two matrix expressions and returns an appropriate
    // matrix expression
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, MatrixExpression<D,E,F> >
    operator*(const MatrixExpression<A,B,C> & a, 
	      const MatrixExpression<D,E,F> & b);

    // operator* takes a matrix expressions and a vector expression, and
    // returns an appropriate vector expression
    template<class A,int B,class C,class D,int E,class F>
    inline VectorExpression < MatrixExpression<A,B,C>, ops::TimesOp, VectorExpression<D,E,F> >
    operator*(const MatrixExpression<A,B,C> & a, 
	      const VectorExpression<D,E,F> & b);

    // operator* takes a scalar and a matrix expression, and returns an
    // appropriate matrix expression
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator* (double b, 
	       const MatrixExpression<A,B,C> & a);

    // operator* takes a matrix expression and a scalar, and returns an
    // appropriate matrix expression
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator* (const MatrixExpression<A,B,C> & a, 
	       double b);

    // operator/ takes a matrix expression and a scalar, and returns an
    // appropriate matrix expression
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator/ (const MatrixExpression<A,B,C> & a, 
	       double b);

    // unary operator- takes a matrix expression and returns an
    // appropriate matrix expression
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::NegativeOp, Base >
    operator- (const MatrixExpression<A,B,C> & a);

    // Transpose function, takes a matrix expression and returns an
    // appropriate matrix expression
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TransposeOp, Base> 
    Transpose(const MatrixExpression<A,B,C> & a);

    // Dyadic function, takes two vector expressions and returns an
    // appropriate matrix expression
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < VectorExpression<A,B,C>, ops::DyadicOp, VectorExpression<D,E,F> > 
    Dyadic(const VectorExpression<A,B,C> & a, 
	   const VectorExpression<D,E,F> & b);

    // Inverse of a matrix
    inline MatrixExpression<> Inverse(const MatrixExpression<> &M);

    // rotation matrix around vector V by an angle V.nrm()
    inline MatrixExpression<> Rodrigues(const Vector &V);

    //                                           
    // member functions of the basic Matrix type 
    //                                           

    // assignment from an expression
    template<class A,int B,class C>
    inline MatrixExpression<>& MatrixExpression<>:: operator=(const MatrixExpression<A,B,C> &e)
    {
      double newvals[3][3] 
	= {{e.eval<0,0>(), e.eval<0,1>(), e.eval<0,2>()},
	   {e.eval<1,0>(), e.eval<1,1>(), e.eval<1,2>()},
	   {e.eval<2,0>(), e.eval<2,1>(), e.eval<2,2>()}};

      xx = newvals[0][0]; xy = newvals[0][1]; xz = newvals[0][2];
      yx = newvals[1][0]; yy = newvals[1][1]; yz = newvals[1][2];
      zx = newvals[2][0]; zy = newvals[2][1]; zz = newvals[2][2];
      return *this;
    }

    // multiply-by assignment from an expression
    inline MatrixExpression<>& MatrixExpression<>:: operator*=(const double d)
    {
      xx *= d;  xy *= d;  xz *= d;
      yx *= d;  yy *= d;  yz *= d;
      zx *= d;  zy *= d;  zz *= d;
      return *this;
    }
    // divide-by assignment from an expression
    inline MatrixExpression<>& MatrixExpression<>:: operator/=(const double d)
    {
      xx /= d;  xy /= d;  xz /= d;
      yx /= d;  yy /= d;  yz /= d;
      zx /= d;  zy /= d;  zz /= d;
      return *this;
    }

    // add-to assignment from an expression
    template<class A,int B,class C>
    inline MatrixExpression<>& MatrixExpression<>:: operator+=(const MatrixExpression<A,B,C> &e)
    { return (*this = *this + e); }
  
    // subtract-from assignment from an expression
    template<class A,int B,class C>
    inline MatrixExpression<>& MatrixExpression<>:: operator-=(const MatrixExpression<A,B,C> &e)
    { return (*this = *this - e); }

    // multiply-by-matrix assignment from an expression
    template<class A,int B,class C>
    inline MatrixExpression<>& MatrixExpression<>::operator*=(const MatrixExpression<A,B,C> & m)
    {
      MatrixExpression<> rhs(*this * m);
      return (*this = rhs);
    }


    // to access the elements of a matrix expression
#define MATPARENTHESIS operator()(const int _i, const int _j) const	\
    {									\
      switch(_i){							\
      case 0: switch(_j){						\
	case 0: return eval<0,0>();					\
	case 1: return eval<0,1>();					\
	case 2: return eval<0,2>();					\
	default: return 0;						\
	}								\
      case 1: switch(_j){						\
	case 0: return eval<1,0>();					\
	case 1: return eval<1,1>();					\
	case 2: return eval<1,2>();					\
	default: return 0;						\
	}								\
      case 2: switch(_j){						\
	case 0: return eval<2,0>();					\
	case 1: return eval<2,1>();					\
	case 2: return eval<2,2>();					\
	default: return 0;						\
	}								\
      default: return 0;						\
      }									\
    }


    // to get the norm of a matrix
#define MATNRM nrm() const					\
    {								\
      double xx=eval<0,0>(), xy=eval<0,1>(), xz=eval<0,2>();	\
      double yx=eval<1,0>(), yy=eval<1,1>(), yz=eval<1,2>();	\
      double zx=eval<2,0>(), zy=eval<2,1>(), zz=eval<2,2>();	\
      double absval;						\
      double biggest=xx<0?-xx:xx;				\
      absval=xy<0?-xy:xy;					\
      if (absval>biggest) biggest=absval;			\
      absval=xz<0?-xz:xz;					\
      if (absval>biggest) biggest=absval;			\
      absval=yx<0?-yx:yx;					\
      if (absval>biggest) biggest=absval;			\
      absval=yy<0?-yy:yy;					\
      if (absval>biggest) biggest=absval;			\
      absval=yz<0?-yz:yz;					\
      if (absval>biggest) biggest=absval;			\
      absval=zx<0?-zx:zx;					\
      if (absval>biggest) biggest=absval;			\
      absval=zy<0?-zy:zy;					\
      if (absval>biggest) biggest=absval;			\
      absval=zz<0?-zz:zz;					\
      if (absval>biggest) biggest=absval;			\
      if (biggest!=0)						\
	return biggest*(double)(sqrt(SQR(xx/biggest)		\
				     +SQR(xy/biggest)		\
				     +SQR(xz/biggest)		\
				     +SQR(yx/biggest)		\
				     +SQR(yy/biggest)		\
				     +SQR(yz/biggest)		\
				     +SQR(zx/biggest)		\
				     +SQR(zy/biggest)		\
				     +SQR(zz/biggest)));	\
      else							\
	return 0;						\
    } 

    // to get the norm squared of a matrix expression
#define MATNRM2 nrm2() const					\
    {								\
      return SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>())	\
        +SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>())	\
        +SQR(eval<0,0>())+SQR(eval<0,0>())+SQR(eval<0,0>());	\
    }

    // to get the trace of a matrix expression
#define MATTR tr() const				\
    {							\
      return eval<0,0>() + eval<1,1>() + eval<2,2>();	\
    }

    // to get the determinant of a matrix expression
#define MATDET det() const						\
    {									\
      double xx=eval<0,0>(), xy=eval<0,1>(), xz=eval<0,2>();		\
      double yx=eval<1,0>(), yy=eval<1,1>(), yz=eval<1,2>();		\
      double zx=eval<2,0>(), zy=eval<2,1>(), zz=eval<2,2>();		\
      return xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx);	\
    }

    // to get the i-th row of a matrix expression
    // note: the class of the matrix expression has to be defined in CLASS
#define MATROW VectorExpression<CLASS,ops::RowOp,Base> row(int i) const \
    {									\
      return VectorExpression<CLASS,ops::RowOp,Base>( *this, i );	\
    }

    // to get the j-th column of a matrix expression
    // note: the class of the matrix expression has to be defined in CLASS
#define MATCOLUMN VectorExpression<CLASS,ops::ColumnOp,Base> column(int j) const \
    {									\
      return VectorExpression<CLASS,ops::ColumnOp,Base>( *this, j );	\
    }


    // Expression template class for '+' between two MatrixExpressions
    //
    // A, C, D, and F are classes, B and E are Operators
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Matrix
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Matrix
    template<class A,int B,class C,class D,int E,class F>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::PlusOp, MatrixExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(), element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS 

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,ops::PlusOp,MatrixExpression<D,E,F> >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Matrix
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Matrix
    template<class A,int B,class C,class D,int E,class F>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::MinusOp, MatrixExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(), and
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,ops::MinusOp,MatrixExpression<D,E,F> >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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


    // Expression template class for '*' between a MatrixExpression and a double
    //
    // A and C are classes, B is an operator
    // if B=ops::NoOp and A=C=Base, the sub-expression is actually a Matrix
    template<class A,int B,class C>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    {
    public:
  
      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,ops::TimesOp,double >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // constructor
      inline 
      MatrixExpression(const MatrixExpression<A,B,C> & _a, 
		       double _b) 
	: a(&_a), b(_b) {}

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
      { return a->eval<I,J>() * b; }

      inline MatrixExpression & operator* (double c)
      {
	b *= c;
	return *this;
      }

      inline MatrixExpression & operator/ (double c)
      {
	b /= c;
	return *this;
      }

    private:
      // store the double and a pointer to the sub-MatrixExpression
      const MatrixExpression<A,B,C>* a;
      double b;

      template<class D,int E,class F> 
      friend MatrixExpression < MatrixExpression<D,E,F>, ops::TimesOp, double >&
      operator* (double b, MatrixExpression < MatrixExpression<D,E,F>, ops::TimesOp, double > &a);

      template<class D,int E,class F> 
      friend MatrixExpression < MatrixExpression<D,E,F>, ops::TimesOp, double >&
      operator/ (double b, MatrixExpression < MatrixExpression<D,E,F>, ops::TimesOp, double > &a);
    };


    // Expression template class for '*' between two MatrixExpressions
    //
    // A, C, D, and F are classes, B and E are Operators
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Matrix
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Matrix
    template<class A,int B,class C,class D,int E,class F>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, MatrixExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression<MatrixExpression<A,B,C>,ops::TimesOp,MatrixExpression<D,E,F> >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the sub-expression is actually a Matrix
    template<class A,int B,class C>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::NegativeOp, Base >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< MatrixExpression<A,B,C>, ops::NegativeOp, Base >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the sub-expression is actually a Matrix
    template<class A,int B,class C>
    class MatrixExpression < MatrixExpression<A,B,C>, ops::TransposeOp, Base >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< MatrixExpression<A,B,C>, ops::TransposeOp, Base >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Vector
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Vector
    template<class A,int B,class C,class D,int E,class F>
    class MatrixExpression < VectorExpression<A,B,C>, ops::DyadicOp, VectorExpression<D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm(),(*this).te(), (*this).det(),
      // element access (*this)(i,j)
      inline double MATNRM2
      inline double MATNRM
      inline double MATDET
      inline double MATTR
      inline double MATPARENTHESIS

      // define row and column access as if they were Vectors
      // CLASS is used in the MATROW and MATCOLUMN macros
#define CLASS MatrixExpression< VectorExpression<A,B,C>, ops::DyadicOp, VectorExpression<D,E,F> >
      inline MATROW
      inline MATCOLUMN
#undef CLASS

      // define what this operation evaluates to  
      template <int I, int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the sub-expression is actually a Matrix
    template<class A,int B,class C>
    class VectorExpression < MatrixExpression<A,B,C>, ops::RowOp, Base >
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to  
      template <const int J> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the sub-expression is actually a Matrix
    template<class A,int B,class C>
    class VectorExpression < MatrixExpression<A,B,C>, ops::ColumnOp, Base>
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to  
      template <const int I> inline double eval() const
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
    // if B=ops::NoOp and A=C=Base, the first sub-expression is actually a Matrix
    // if E=ops::NoOp and D=F=Base, the second sub-expression is actually a Vector
    template<class A,int B,class C,class D,int E,class F>
    class VectorExpression < MatrixExpression<A,B,C>, ops::TimesOp, VectorExpression <D,E,F> >
    {
    public:

      // define (*this).nrm2(), (*this).nrm() and element access (*this)(i) 
      inline double VECNRM2
      inline double VECNRM
      inline double VECPARENTHESES

      // define what this operation evaluates to  
      template <int I> inline double eval() const
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

    // Matrix + Matrix
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::PlusOp, MatrixExpression<D,E,F> >
    operator+ (const MatrixExpression<A,B,C> & a, 
	       const MatrixExpression<D,E,F> & b)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::PlusOp, MatrixExpression<D,E,F> >(a,b);
    }

    // Matrix - Matrix
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::MinusOp, MatrixExpression<D,E,F> >
    operator- (const MatrixExpression<A,B,C> & a, 
	       const MatrixExpression<D,E,F> & b)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::MinusOp, MatrixExpression<D,E,F> >(a, b);
    }

    // double * Matrix
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator* (double b, 
	       const MatrixExpression<A,B,C> & a)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >(a, b);
    }

    // Matrix * double
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator* (const MatrixExpression<A,B,C> & a, 
	       double b)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >(a, b);
    }

    // Matrix / double
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >
    operator/ (const MatrixExpression<A,B,C> & a, 
	       double b)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >(a, 1.0/b);
    }

    // Matrix * double * double
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >&
    operator* (double b, MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double > &a)
    { 
      a.b *= b; 
      return a; 
    }

    // Matrix * double / double
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double >&
    operator/ (double b, MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, double > &a)
    { 
      a.b /= b; 
      return a; 
    }

    // Matrix * Matrix
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, MatrixExpression<D,E,F> >
    operator*(const MatrixExpression<A,B,C> & a, 
	      const MatrixExpression<D,E,F> & b)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::TimesOp, MatrixExpression<D,E,F> > (a, b);
    }

    // -Matrix 
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::NegativeOp, Base >
    operator- (const MatrixExpression<A,B,C> & a)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::NegativeOp, Base > (a);
    }

    // Matrix * Vector
    template<class A,int B,class C,class D,int E,class F>
    inline VectorExpression < MatrixExpression<A,B,C>, ops::TimesOp, VectorExpression<D,E,F> >
    operator*(const MatrixExpression<A,B,C> & a, 
	      const VectorExpression<D,E,F> & b)
    { 
      return VectorExpression < MatrixExpression<A,B,C>, ops::TimesOp, VectorExpression<D,E,F> > (a, b);
    }

    //
    // Implementation of the row and column member functions of MatrixExpression
    //

    inline VectorExpression<MatrixExpression<>,ops::RowOp,Base> 
    MatrixExpression<>::row(int i)
    {
      return VectorExpression<MatrixExpression<>,ops::RowOp,Base>(*this,i);
    }

    inline VectorExpression<MatrixExpression<>,ops::ColumnOp,Base> 
    MatrixExpression<>::column(int j)
    {
      return VectorExpression<MatrixExpression<>,ops::ColumnOp,Base>(*this,j);
    }

    // Can now implement Gramm-Schmidt orthogonalization of the rows of the matrix
    inline void Matrix::reorthogonalize()
    {
      double z;
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

    // the inverse of a matrix
    inline MatrixExpression<> Inverse(const MatrixExpression<> &M) 
    {
      double d = 1/M.det();
      return MatrixExpression<>(
				(M.yy*M.zz-M.yz*M.zy)*d, -(M.xy*M.zz-M.xz*M.zy)*d,  (M.xy*M.yz-M.xz*M.yy)*d,
				-(M.yx*M.zz-M.yz*M.zx)*d,  (M.xx*M.zz-M.xz*M.zx)*d, -(M.xx*M.yz-M.xz*M.yx)*d,
				(M.yx*M.zy-M.yy*M.zx)*d, -(M.xx*M.zy-M.xy*M.zx)*d,  (M.xx*M.yy-M.xy*M.yx)*d);
    }

    // rotation matrix built using the Rodrigues formula
    inline MatrixExpression<> Rodrigues(const Vector &V) 
    {
      double theta = V.nrm();
      if (theta != 0) {
	double s(std::sin(theta)), c(std::cos(theta));
	double inrm = 1/theta;
	double wx = V(0)*inrm;
	double wy = V(1)*inrm;
	double wz = V(2)*inrm;
	double oneminusc = 1-c;
	double wxwy1mc = wx*wy*oneminusc;
	double wxwz1mc = wx*wz*oneminusc;
	double wywz1mc = wy*wz*oneminusc;
	double wxs = wx*s;
	double wys = wy*s;
	double wzs = wz*s;
	return MatrixExpression<>(c+wx*wx*oneminusc,wxwy1mc-wzs,      wxwz1mc+wys,
				  wxwy1mc+wzs,      c+wy*wy*oneminusc,wywz1mc-wxs,
				  wxwz1mc-wys,     wywz1mc+wxs, c+wz*wz*oneminusc);
      } else {
	return MatrixExpression<>(1,0,0,0,1,0,0,0,1);
      }
    }

    // transpose of a matrix
    template<class A,int B,class C>
    inline MatrixExpression < MatrixExpression<A,B,C>, ops::TransposeOp, Base >
    Transpose (const MatrixExpression<A,B,C> & a)
    { 
      return MatrixExpression < MatrixExpression<A,B,C>, ops::TransposeOp, Base > (a);
    }

    // dyadic product of two vectors
    template<class A,int B,class C,class D,int E,class F>
    inline MatrixExpression < VectorExpression<A,B,C>, ops::DyadicOp, VectorExpression<D,E,F> >
    Dyadic (const VectorExpression<A,B,C> & a, 
	    const VectorExpression<D,E,F> & b)
    { 
      return MatrixExpression < VectorExpression<A,B,C>, ops::DyadicOp, VectorExpression<D,E,F> > (a,b);
    }

#undef MATNRM
#undef MATNRM2
#undef MATPARENTHESIS
#undef MATTR
#undef MATDET
#undef MATROW
#undef MATCOLUMN

    // vectors
    template<class A, int B, class C>
    inline magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
					      const MatrixExpression<A,B,C> & t )
    {
      char name[2] = "x";

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	{
	  name[0] = 'x'+iDim;
	  XML << magnet::xml::tag(name);
      
	  for (size_t jDim = 0; jDim < NDIM; ++jDim)
	    {
	      char name2[2] = "x";
	      name2[0] = 'x'+jDim;
	      XML << magnet::xml::attr(name2) << t(iDim,jDim);
	    }

	  XML << magnet::xml::endtag(name);
	}
  
      return XML;
    }

    inline
    MatrixExpression<>& 
    operator<<(MatrixExpression<>& data, const magnet::xml::Node& XML)
    {
      char name[2] = "x";
  
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	{
	  name[0] = 'x'+iDim;

	  for (size_t jDim = 0; jDim < NDIM; ++jDim)
	    {
	      char name2[2] = "x";
	      name2[0] = 'x'+jDim;

	      try {
		data(iDim,jDim) = XML.getNode(name).getAttribute(name2).as<double>();
	      }
	      catch (boost::bad_lexical_cast &)
		{
		  M_throw() << "Failed a lexical cast";
		}
	    }
	}

      return data;
    }
  }
}

namespace coil { typedef ::magnet::math::Matrix Matrix; }
namespace dynamo { typedef ::magnet::math::Matrix Matrix; }

#define BOOST_TEST_MODULE Vector_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/vector.hpp>
#include <magnet/math/matrix.hpp>
#include <cmath>

std::mt19937 RNG;
std::normal_distribution<double> normal_dist(0, 1);
std::uniform_real_distribution<double> angle_dist(0, std::atan(1)*4);
using namespace magnet::math;

Vector random_unit_vec() {
  Vector vec{normal_dist(RNG), normal_dist(RNG), normal_dist(RNG)};
  return vec / vec.nrm();
}

using namespace magnet::math;
BOOST_AUTO_TEST_CASE( vector_initializer_list )
{
  Vector A{1,2,3};
  BOOST_CHECK_EQUAL(A[0], 1);
  BOOST_CHECK_EQUAL(A[1], 2);
  BOOST_CHECK_EQUAL(A[2], 3);
}

BOOST_AUTO_TEST_CASE( vector_assignment )
{
  Vector A{1,2,3};
  Vector B;
  B = A;
  BOOST_CHECK_EQUAL(B[0], 1);
  BOOST_CHECK_EQUAL(B[1], 2);
  BOOST_CHECK_EQUAL(B[2], 3);
}

BOOST_AUTO_TEST_CASE( vector_comparison )
{
  Vector A{1,2,3};
  Vector B{4,5,6};
  BOOST_CHECK(A!=B);
  BOOST_CHECK(B!=A);
  BOOST_CHECK(A==A);
  BOOST_CHECK(B==B);
}

BOOST_AUTO_TEST_CASE( vector_addition )
{
  Vector A{1,2,3};
  Vector B{4,5,6};

  Vector C = A+B;
  BOOST_CHECK_EQUAL(C[0], 5);
  BOOST_CHECK_EQUAL(C[1], 7);
  BOOST_CHECK_EQUAL(C[2], 9);

  Vector D=A;
  D+=B;
  BOOST_CHECK_EQUAL(D[0], 5);
  BOOST_CHECK_EQUAL(D[1], 7);
  BOOST_CHECK_EQUAL(D[2], 9);
}

BOOST_AUTO_TEST_CASE( vector_subtraction )
{
  Vector A{1,2,3};
  Vector B{4,5,6};

  Vector C = A-B;
  BOOST_CHECK_EQUAL(C[0], -3);
  BOOST_CHECK_EQUAL(C[1], -3);
  BOOST_CHECK_EQUAL(C[2], -3);

  Vector D=A;
  D-=B;
  BOOST_CHECK_EQUAL(D[0], -3);
  BOOST_CHECK_EQUAL(D[1], -3);
  BOOST_CHECK_EQUAL(D[2], -3);
}

BOOST_AUTO_TEST_CASE( vector_scalar_prod )
{
  Vector A{1,2,3};
  Vector B{4,5,6};

  BOOST_CHECK_EQUAL((A*B), 32);
  BOOST_CHECK_EQUAL((A|B), 32);
  BOOST_CHECK_EQUAL((Vector{1,1,0} * Vector{0,0,1}), 0);
}

BOOST_AUTO_TEST_CASE( vector_cross_prod )
{
  Vector A{1,2,3};
  Vector B{3,2,1};
  Vector C=A^B;
  BOOST_CHECK(C == (Vector{-4,8,-4}));
}

BOOST_AUTO_TEST_CASE( vector_unary_negative )
{
  Vector A{1,2,3};
  Vector B(-A);

  BOOST_CHECK_EQUAL(B[0], -1);
  BOOST_CHECK_EQUAL(B[1], -2);
  BOOST_CHECK_EQUAL(B[2], -3);
}

BOOST_AUTO_TEST_CASE( vector_float_mult )
{
  Vector A{1,2,3};
  Vector C = A*2.0;
  BOOST_CHECK_EQUAL(C[0], 2);
  BOOST_CHECK_EQUAL(C[1], 4);
  BOOST_CHECK_EQUAL(C[2], 6);
  C = 2.0 * A;
  BOOST_CHECK_EQUAL(C[0], 2);
  BOOST_CHECK_EQUAL(C[1], 4);
  BOOST_CHECK_EQUAL(C[2], 6);
}

BOOST_AUTO_TEST_CASE( vector_norm )
{
  Vector B{1,1,1};
  BOOST_CHECK_CLOSE(B.nrm2(), 3, 1e-10);
  BOOST_CHECK_CLOSE(B.nrm(), std::sqrt(3), 1e-10);
  B.normalise();
  BOOST_CHECK_EQUAL(B[0], 1.0/std::sqrt(3.0));
  BOOST_CHECK_EQUAL(B[1], 1.0/std::sqrt(3.0));
  BOOST_CHECK_EQUAL(B[2], 1.0/std::sqrt(3.0));

  B = Vector{-1,0,0};
  BOOST_CHECK_CLOSE(B.nrm2(), 1, 1e-10);
  BOOST_CHECK_CLOSE(B.nrm(), 1, 1e-10);
  B.normalise();
  BOOST_CHECK_CLOSE(B[0], -1.0, 1e-10);
  BOOST_CHECK_CLOSE(B[1], 0, 1e-10);
  BOOST_CHECK_CLOSE(B[2], 0, 1e-10);

  B = Vector{0,0,0};
  BOOST_CHECK_EQUAL(B.nrm2(), 0);
  BOOST_CHECK_EQUAL(B.nrm(), 0);
  B.normalise();
  BOOST_CHECK_EQUAL(B[0], 0);
  BOOST_CHECK_EQUAL(B[1], 0);
  BOOST_CHECK_EQUAL(B[2], 0);
}

BOOST_AUTO_TEST_CASE( matrix_identity )
{
  Matrix B{1,0,0,0,1,0,0,0,1};
  BOOST_CHECK(B == Matrix::identity());
}

BOOST_AUTO_TEST_CASE( matrix_comparison )
{
  Matrix B{1,2,3,4,1,6,7,8,1};
  BOOST_CHECK(Matrix::identity() == Matrix::identity());
  BOOST_CHECK(B != Matrix::identity());
  BOOST_CHECK(B == B);
}

BOOST_AUTO_TEST_CASE( matrix_matrix_multiplication )
{
  BOOST_CHECK(Matrix::identity() * Matrix::identity() == Matrix::identity());
  
  Matrix A{1,2,3,4,5,6,7,8,9};
  Matrix B = A * A;
  Matrix result{30,36,42,66,81,96,102,126,150};
  BOOST_CHECK(B == result);

  B=Matrix{3,2,1,4,5,6,9,8,7};
  Matrix C = A * B;
  result=Matrix{38,36,34,86,81,76,134,126,118};
  BOOST_CHECK(C == result);  
}

BOOST_AUTO_TEST_CASE( matrix_matrix_multiplication_4D)
{
  typedef NMatrix<double, 4> Matrix;

  BOOST_CHECK(Matrix::identity() * Matrix::identity() == Matrix::identity());

  const Matrix A{1,0,1,-2,0,1,0,2,2,0,1,0,-1,1,0,1};
  Matrix B = Matrix::identity() * A;
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(B(i), A(i), 0.0000001);
}

BOOST_AUTO_TEST_CASE( matrix_vector_multiplication )
{
  RNG.seed(1);
  const Vector vec = random_unit_vec();
  Vector result = Matrix::identity() * vec;
  BOOST_CHECK(vec == result);
  result = (-Matrix::identity()) * vec;
  BOOST_CHECK(-vec == result);

  Matrix A{1,2,3,4,5,6,7,8,9};
  Vector b{2,3,4};
  Vector c = A * b;
  BOOST_CHECK((c == Vector{20,47,74}));
}

BOOST_AUTO_TEST_CASE( matrix_scalar_multiplication )
{
  const Matrix A{1,2,3,4,5,6,7,8,9};
  Matrix B= A * 2;
  Matrix C= A * 2.0;
  BOOST_CHECK(B == C);
  B= 2 * A;
  C= 2.0 * A;
  BOOST_CHECK(B == C);

  B=A;
  B*=2;
  C=A;
  C*=2.0;
  BOOST_CHECK(B == C);
  for (size_t i(0); i < 9; ++i)
    BOOST_CHECK(2.0 * A(i) == C(i));
}

BOOST_AUTO_TEST_CASE( matrix_dyadic )
{
  const Vector A{1,2,3};
  const Vector B{4,5,6};
  const Matrix C = Dyadic(A,B);
  const Matrix result{4,5,6,8,10,12,12,15,18};
  BOOST_CHECK(C == result);
}

BOOST_AUTO_TEST_CASE( matrix_determinant_2D )
{
  const NMatrix<double, 2> A{1,2,3,4};
  BOOST_CHECK_CLOSE(determinant(A), -2, 1e-10);
}

BOOST_AUTO_TEST_CASE( matrix_determinant_3D )
{
  const Matrix A{1,0,3,4,5,6,9,8,7};
  BOOST_CHECK_CLOSE(determinant(A), -52, 1e-10);
}

BOOST_AUTO_TEST_CASE( matrix_determinant_4D )
{
  typedef NMatrix<double, 4> Matrix;
  const Matrix A{1,0,1,-2,0,1,0,2,2,0,1,0,-1,1,0,1};
  BOOST_CHECK_CLOSE(determinant(A), -1, 1e-10);
}

BOOST_AUTO_TEST_CASE( matrix_inverse_2D )
{
  typedef NMatrix<double, 2> Matrix;
  const Matrix A{1,2,3,4};
  const Matrix B = inverse(A);
  const Matrix result{-2, 1, 3.0/2, -0.5};

  Matrix I = result * A;
  for (size_t i(0); i < 2*2; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  I = A * result;
  for (size_t i(0); i < 2*2; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  for (size_t i(0); i < 2*2; ++i)
    BOOST_CHECK_CLOSE(result(i), B(i), 0.0000001);

  I = A * B;
  for (size_t i(0); i < 2*2; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  I = B * A;
  for (size_t i(0); i < 2*2; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
}


BOOST_AUTO_TEST_CASE( matrix_inverse_3D )
{
  const Matrix A{1,1,3,0,1,3,1,0,1};
  const Matrix B = inverse(A);
  const Matrix result{1,-1,0,3,-2,-3,-1,1,1};

  Matrix I = result * A;
  for (size_t i(0); i < 3*3; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  I = A * result;
  for (size_t i(0); i < 3*3; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  for (size_t i(0); i < 3*3; ++i)
    BOOST_CHECK_CLOSE(result(i), B(i), 0.0000001);

  I = A * B;
  for (size_t i(0); i < 3*3; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);

  I = B * A;
  for (size_t i(0); i < 3*3; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
}

BOOST_AUTO_TEST_CASE( matrix_inverse_4D )
{
  typedef NMatrix<double, 4> Matrix;
  const Matrix A{1,0,1,-2,0,1,0,2,2,0,1,0,-1,1,0,1};
  const Matrix B = inverse(A);
  const Matrix result{1,2,-1,-2,2,3,-2,-2,-2,-4,3,4,-1,-1,1,1};

  Matrix I = result * A;
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
  
  I = A * result;
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
  
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(result(i), B(i), 0.0000001);
  
  I = A * B;
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
  
  I = B * A;
  for (size_t i(0); i < 4*4; ++i)
    BOOST_CHECK_CLOSE(Matrix::identity()(i), I(i), 0.0000001);
}

#include <magnet/math/polynomial.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/math/trigsymbols.hpp>

const size_t testcount = 100;
const double errlvl = 1e-10;

BOOST_AUTO_TEST_CASE( Vector_symbolic )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  //A tough test is to implement the Rodriugues formula symbolically.
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      double angle = angle_dist(RNG);
      Vector axis = random_unit_vec();
      Vector start = random_unit_vec();
      Vector end = Rodrigues(axis * angle) * start;
      
      Vector r = axis * (axis * start);
      auto f = (start - r) * cos(x) + (axis ^ start) * sin(x) + r;
      Vector err = end - eval(f, angle);
      
      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

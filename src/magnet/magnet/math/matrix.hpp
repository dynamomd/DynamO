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
#include <magnet/math/detail/eigenval.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
namespace math {
/*! \brief N-dimensional matrix type.

  \tparam N The dimensionality of the NMatrix.

  \tparam T The type of the components of the
  NMatrix.
 */
template <class T = double, size_t N = 3>
class NMatrix : public std::array<T, N * N> {
  typedef std::array<T, N * N> Base;

public:
  inline NMatrix() { Base::fill(T()); }

  NMatrix(std::initializer_list<T> _list) {
    if (_list.size() > N * N)
      throw std::length_error("initializer list too long");

    size_t i = 0;
    auto it = _list.begin();
    for (; it != _list.end(); ++i, ++it)
      operator()(i) = *it;

    for (; i < N * N; ++i)
      operator()(i) = 0.0;
  }

  inline T tr() const {
    T sum(0);
    for (size_t i(0); i < N; ++i)
      sum += operator()(i, i);
    return sum;
  }

  static inline NMatrix<T, N> identity() {
    NMatrix retval;
    for (size_t i(0); i < N; ++i)
      retval(i, i) = T(1);
    return retval;
  }

  NMatrix<T, N> transpose() const {
    NMatrix<T, N> retval;
    for (size_t i(0); i < N; ++i)
      for (size_t j(0); j < N; ++j)
        retval(i, j) = operator()(j, i);
    return retval;
  }

  inline T &operator()(size_t i, size_t j) {
    return *(Base::begin() + N * i + j);
  }
  inline const T &operator()(size_t i, size_t j) const {
    return *(Base::begin() + N * i + j);
  }

  inline T &operator()(size_t i) { return *(Base::begin() + i); }
  inline const T &operator()(size_t i) const { return *(Base::begin() + i); }

  inline bool operator==(const NMatrix &omat) const {
    for (size_t i(0); i < N * N; ++i)
      if (operator()(i) != omat(i))
        return false;
    return true;
  }

  NMatrix<T, N> &operator*=(const T d) {
    for (size_t i(0); i < N * N; ++i)
      Base::operator[](i) *= d;
    return *this;
  }

  NMatrix<T, N> &operator/=(const T d) {
    for (size_t i(0); i < N * N; ++i)
      Base::operator[](i) /= d;
    return *this;
  }

  NMatrix<T, N> &operator+=(const NMatrix<T, N> M) {
    for (size_t i(0); i < N * N; ++i)
      Base::operator[](i) += M(i);
    return *this;
  }

  NMatrix<T, N> &operator-=(const NMatrix<T, N> M) {
    for (size_t i(0); i < N * N; ++i)
      Base::operator[](i) -= M(i);
    return *this;
  }

  NMatrix<T, N> &operator*=(const NMatrix<T, N> M) { return *this = *this * M; }

  NMatrix<T, N> operator-() const {
    NMatrix<T, N> retval;
    for (size_t i(0); i < N * N; ++i)
      retval(i) = -operator()(i);
    return retval;
  }

  std::string toString() const {
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<double>::digits10 + 2) << "<";
    for (size_t i(0); i < N * N - 1; ++i)
      os << Base::operator[](i) << ",";
    os << Base::operator[](N * N - 1) << ">";
    return os.str();
  }
};

/*! \relates NMatrix
  \name Determinant
  Methods for calculating the determinant of low-order NMatrix types.
  \{
 */
template <class T, size_t N> T determinant(const NMatrix<T, N> &M) {
  static_assert(!N, "Determinant not implemented for this dimension of Matrix");
}

/*! \brief Determinant of a 1x1 NMatrix. */
template <class T> T determinant(const NMatrix<T, 1> &M) { return M(0, 0); }

/*! \brief Determinant of a 2x2 NMatrix. */
template <class T> T determinant(const NMatrix<T, 2> &M) {
  return M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
}

/*! \brief Determinant of a 3x3 NMatrix. */
template <class T> T determinant(const NMatrix<T, 3> &M) {
  return M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) +
         M(0, 1) * (M(1, 2) * M(2, 0) - M(1, 0) * M(2, 2)) +
         M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
}

/*! \brief Determinant of a 4x4 NMatrix. */
template <class T> T determinant(const NMatrix<T, 4> &M) {
  return M(0, 3) * M(1, 2) * M(2, 1) * M(3, 0) -
         M(0, 2) * M(1, 3) * M(2, 1) * M(3, 0) -
         M(0, 3) * M(1, 1) * M(2, 2) * M(3, 0) +
         M(0, 1) * M(1, 3) * M(2, 2) * M(3, 0) +
         M(0, 2) * M(1, 1) * M(2, 3) * M(3, 0) -
         M(0, 1) * M(1, 2) * M(2, 3) * M(3, 0) -
         M(0, 3) * M(1, 2) * M(2, 0) * M(3, 1) +
         M(0, 2) * M(1, 3) * M(2, 0) * M(3, 1) +
         M(0, 3) * M(1, 0) * M(2, 2) * M(3, 1) -
         M(0, 0) * M(1, 3) * M(2, 2) * M(3, 1) -
         M(0, 2) * M(1, 0) * M(2, 3) * M(3, 1) +
         M(0, 0) * M(1, 2) * M(2, 3) * M(3, 1) +
         M(0, 3) * M(1, 1) * M(2, 0) * M(3, 2) -
         M(0, 1) * M(1, 3) * M(2, 0) * M(3, 2) -
         M(0, 3) * M(1, 0) * M(2, 1) * M(3, 2) +
         M(0, 0) * M(1, 3) * M(2, 1) * M(3, 2) +
         M(0, 1) * M(1, 0) * M(2, 3) * M(3, 2) -
         M(0, 0) * M(1, 1) * M(2, 3) * M(3, 2) -
         M(0, 2) * M(1, 1) * M(2, 0) * M(3, 3) +
         M(0, 1) * M(1, 2) * M(2, 0) * M(3, 3) +
         M(0, 2) * M(1, 0) * M(2, 1) * M(3, 3) -
         M(0, 0) * M(1, 2) * M(2, 1) * M(3, 3) -
         M(0, 1) * M(1, 0) * M(2, 2) * M(3, 3) +
         M(0, 0) * M(1, 1) * M(2, 2) * M(3, 3);
}
/*! \} */

/*! \relates NMatrix
  \name Adjoint
  Methods for calculating the adjoint of low-order NMatrix types.
  \{
 */
template <class T, size_t N> NMatrix<T, N> adjoint(const NMatrix<T, N> &M) {
  static_assert(!N, "Adjoint not implemented for this dimension of Matrix");
}

/*! \brief Adjoint of a 1x1 NMatrix. */
template <class T> NMatrix<T, 1> adjoint(const NMatrix<T, 1> &M) {
  return NMatrix<T, 1>{1};
}

/*! \brief Adjoint of a 2x2 NMatrix. */
template <class T> NMatrix<T, 2> adjoint(const NMatrix<T, 2> &M) {
  return NMatrix<T, 2>{M(1, 1), -M(0, 1), -M(1, 0), M(0, 0)};
}

/*! \brief Adjoint of a 3x3 NMatrix. */
template <class T> NMatrix<T, 3> adjoint(const NMatrix<T, 3> &M) {
  return NMatrix<T, 3>{(M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)),
                       -(M(0, 1) * M(2, 2) - M(0, 2) * M(2, 1)),
                       (M(0, 1) * M(1, 2) - M(0, 2) * M(1, 1)),
                       -(M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)),
                       (M(0, 0) * M(2, 2) - M(0, 2) * M(2, 0)),
                       -(M(0, 0) * M(1, 2) - M(0, 2) * M(1, 0)),
                       (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0)),
                       -(M(0, 0) * M(2, 1) - M(0, 1) * M(2, 0)),
                       (M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0))};
}

/*! \brief Adjoint of a 4x4 NMatrix. */
template <class T> NMatrix<T, 4> adjoint(const NMatrix<T, 4> &M) {
  NMatrix<T, 4> R;
  R(0, 0) = M(1, 2) * M(2, 3) * M(3, 1) - M(1, 3) * M(2, 2) * M(3, 1) +
            M(1, 3) * M(2, 1) * M(3, 2) - M(1, 1) * M(2, 3) * M(3, 2) -
            M(1, 2) * M(2, 1) * M(3, 3) + M(1, 1) * M(2, 2) * M(3, 3);
  R(0, 1) = M(0, 3) * M(2, 2) * M(3, 1) - M(0, 2) * M(2, 3) * M(3, 1) -
            M(0, 3) * M(2, 1) * M(3, 2) + M(0, 1) * M(2, 3) * M(3, 2) +
            M(0, 2) * M(2, 1) * M(3, 3) - M(0, 1) * M(2, 2) * M(3, 3);
  R(0, 2) = M(0, 2) * M(1, 3) * M(3, 1) - M(0, 3) * M(1, 2) * M(3, 1) +
            M(0, 3) * M(1, 1) * M(3, 2) - M(0, 1) * M(1, 3) * M(3, 2) -
            M(0, 2) * M(1, 1) * M(3, 3) + M(0, 1) * M(1, 2) * M(3, 3);
  R(0, 3) = M(0, 3) * M(1, 2) * M(2, 1) - M(0, 2) * M(1, 3) * M(2, 1) -
            M(0, 3) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 3) * M(2, 2) +
            M(0, 2) * M(1, 1) * M(2, 3) - M(0, 1) * M(1, 2) * M(2, 3);
  R(1, 0) = M(1, 3) * M(2, 2) * M(3, 0) - M(1, 2) * M(2, 3) * M(3, 0) -
            M(1, 3) * M(2, 0) * M(3, 2) + M(1, 0) * M(2, 3) * M(3, 2) +
            M(1, 2) * M(2, 0) * M(3, 3) - M(1, 0) * M(2, 2) * M(3, 3);
  R(1, 1) = M(0, 2) * M(2, 3) * M(3, 0) - M(0, 3) * M(2, 2) * M(3, 0) +
            M(0, 3) * M(2, 0) * M(3, 2) - M(0, 0) * M(2, 3) * M(3, 2) -
            M(0, 2) * M(2, 0) * M(3, 3) + M(0, 0) * M(2, 2) * M(3, 3);
  R(1, 2) = M(0, 3) * M(1, 2) * M(3, 0) - M(0, 2) * M(1, 3) * M(3, 0) -
            M(0, 3) * M(1, 0) * M(3, 2) + M(0, 0) * M(1, 3) * M(3, 2) +
            M(0, 2) * M(1, 0) * M(3, 3) - M(0, 0) * M(1, 2) * M(3, 3);
  R(1, 3) = M(0, 2) * M(1, 3) * M(2, 0) - M(0, 3) * M(1, 2) * M(2, 0) +
            M(0, 3) * M(1, 0) * M(2, 2) - M(0, 0) * M(1, 3) * M(2, 2) -
            M(0, 2) * M(1, 0) * M(2, 3) + M(0, 0) * M(1, 2) * M(2, 3);
  R(2, 0) = M(1, 1) * M(2, 3) * M(3, 0) - M(1, 3) * M(2, 1) * M(3, 0) +
            M(1, 3) * M(2, 0) * M(3, 1) - M(1, 0) * M(2, 3) * M(3, 1) -
            M(1, 1) * M(2, 0) * M(3, 3) + M(1, 0) * M(2, 1) * M(3, 3);
  R(2, 1) = M(0, 3) * M(2, 1) * M(3, 0) - M(0, 1) * M(2, 3) * M(3, 0) -
            M(0, 3) * M(2, 0) * M(3, 1) + M(0, 0) * M(2, 3) * M(3, 1) +
            M(0, 1) * M(2, 0) * M(3, 3) - M(0, 0) * M(2, 1) * M(3, 3);
  R(2, 2) = M(0, 1) * M(1, 3) * M(3, 0) - M(0, 3) * M(1, 1) * M(3, 0) +
            M(0, 3) * M(1, 0) * M(3, 1) - M(0, 0) * M(1, 3) * M(3, 1) -
            M(0, 1) * M(1, 0) * M(3, 3) + M(0, 0) * M(1, 1) * M(3, 3);
  R(2, 3) = M(0, 3) * M(1, 1) * M(2, 0) - M(0, 1) * M(1, 3) * M(2, 0) -
            M(0, 3) * M(1, 0) * M(2, 1) + M(0, 0) * M(1, 3) * M(2, 1) +
            M(0, 1) * M(1, 0) * M(2, 3) - M(0, 0) * M(1, 1) * M(2, 3);
  R(3, 0) = M(1, 2) * M(2, 1) * M(3, 0) - M(1, 1) * M(2, 2) * M(3, 0) -
            M(1, 2) * M(2, 0) * M(3, 1) + M(1, 0) * M(2, 2) * M(3, 1) +
            M(1, 1) * M(2, 0) * M(3, 2) - M(1, 0) * M(2, 1) * M(3, 2);
  R(3, 1) = M(0, 1) * M(2, 2) * M(3, 0) - M(0, 2) * M(2, 1) * M(3, 0) +
            M(0, 2) * M(2, 0) * M(3, 1) - M(0, 0) * M(2, 2) * M(3, 1) -
            M(0, 1) * M(2, 0) * M(3, 2) + M(0, 0) * M(2, 1) * M(3, 2);
  R(3, 2) = M(0, 2) * M(1, 1) * M(3, 0) - M(0, 1) * M(1, 2) * M(3, 0) -
            M(0, 2) * M(1, 0) * M(3, 1) + M(0, 0) * M(1, 2) * M(3, 1) +
            M(0, 1) * M(1, 0) * M(3, 2) - M(0, 0) * M(1, 1) * M(3, 2);
  R(3, 3) = M(0, 1) * M(1, 2) * M(2, 0) - M(0, 2) * M(1, 1) * M(2, 0) +
            M(0, 2) * M(1, 0) * M(2, 1) - M(0, 0) * M(1, 2) * M(2, 1) -
            M(0, 1) * M(1, 0) * M(2, 2) + M(0, 0) * M(1, 1) * M(2, 2);
  return R;
}

/*! \} */

/*! \relates NMatrix
  \brief Calculation of the inverse of a NMatrix.

  This is only possible if a suitable \ref adjoint and \ref
  determinant are defined.
 */
template <class T, size_t N> NMatrix<T, N> inverse(const NMatrix<T, N> &M) {
  T det = determinant(M);
  if (det == 0)
    return NMatrix<T, N>::identity();
  return adjoint(M) * (1 / det);
}

/*! \relates NMatrix
  \brief Addition of two NMatrix types.
 */
template <class T1, class T2, size_t N>
NMatrix<decltype(T1() + T2()), N> operator+(const NMatrix<T1, N> &A,
                                            const NMatrix<T2, N> &B) {
  NMatrix<decltype(T1() + T2()), N> retval;
  for (size_t i(0); i < N * N; ++i)
    retval(i) = A(i) + B(i);
  return retval;
}

/*! \relates NMatrix
  \brief Subtraction of two NMatrix types.
 */
template <class T1, class T2, size_t N>
NMatrix<decltype(T1() - T2()), N> operator-(const NMatrix<T1, N> &A,
                                            const NMatrix<T2, N> &B) {
  NMatrix<decltype(T1() - T2()), N> retval;
  for (size_t i(0); i < N * N; ++i)
    retval(i) = A(i) - B(i);
  return retval;
}

/*! \relates NMatrix
  \brief Multiplication of a NMatrix type and a Scalar.
 */
template <class T, class U, size_t N>
NMatrix<decltype(T() * U()), N> operator*(const NMatrix<T, N> &A, const U &B) {
  NMatrix<decltype(T() * U()), N> retval;
  for (size_t i(0); i < N * N; ++i)
    retval(i) = A(i) * B;
  return retval;
}

/*! \relates NMatrix
  \brief Multiplication of a Scalar and an NMatrix type.
 */
template <class T, class U, size_t N>
NMatrix<decltype(U() * T()), N> operator*(const U &B, const NMatrix<T, N> &A) {
  NMatrix<decltype(U() * T()), N> retval;
  for (size_t i(0); i < N * N; ++i)
    retval(i) = B * A(i);
  return retval;
}

/*! \relates NMatrix
  \brief Dot product between two NMatrix types.
 */
template <class T1, class T2, size_t N>
NMatrix<decltype(T1() * T2()), N> operator*(const NMatrix<T1, N> &A,
                                            const NMatrix<T2, N> &B) {
  NMatrix<decltype(T1() * T2()), N> retval;
  for (size_t i(0); i < N; ++i)
    for (size_t j(0); j < N; ++j) {
      decltype(T1() * T2()) sum(0);
      for (size_t k(0); k < N; ++k)
        sum += A(i, k) * B(k, j);
      retval(i, j) = sum;
    }
  return retval;
}

/*! \relates NMatrix
  \brief Dot product between a NMatrix and a NVector.
 */
template <class T1, class T2, size_t N>
NVector<decltype(T1() * T2()), N> operator*(const NMatrix<T1, N> &A,
                                            const NVector<T2, N> &B) {
  NVector<decltype(T1() * T2()), N> retval;
  for (size_t i(0); i < N; ++i) {
    decltype(T1() * T2()) sum(0);
    for (size_t j(0); j < N; ++j)
      sum += A(i, j) * B[j];
    retval[i] = sum;
  }
  return retval;
}

/*! \relates NMatrix
  \brief Dot product between a NVector and a NMatrix.
 */
template <class T1, class T2, size_t N>
NVector<decltype(T1() * T2()), N> operator*(const NVector<T1, N> &A,
                                            const NMatrix<T2, N> &B) {
  NVector<decltype(T1() * T2()), N> retval;
  for (size_t i(0); i < N; ++i) {
    decltype(T1() * T2()) sum(0);
    for (size_t j(0); j < N; ++j)
      sum += A(j, i) * B[j];
    retval[i] = sum;
  }
  return retval;
}

/*! \relates NMatrix
  \brief Division of a NMatrix by a scalar.
 */
template <class T, class U, size_t N>
NMatrix<T, N> operator/(const NMatrix<T, N> &A, const U &B) {
  NMatrix<T, N> retval;
  for (size_t i(0); i < N * N; ++i)
    retval(i) = A(i) / B;
  return retval;
}

/*! \relates NMatrix
  \brief Dyadic product of two NVector types to form a NMatrix.
 */
template <class T, size_t N>
NMatrix<T, N> Dyadic(const NVector<T, N> &A, const NVector<T, N> &B) {
  NMatrix<T, N> retval;
  for (size_t i(0); i < N; ++i)
    for (size_t j(0); j < N; ++j) {
      retval(i, j) = A[i] * B[j];
    }
  return retval;
}

/*! \relates NMatrix
  \brief Create a NMatrix which is equivalent to a cross product
  with the supplied NVector.
 */
template <class T> inline NMatrix<T, 3> cross(const NVector<T, 3> &v) {
  return NMatrix<T, 3>{0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0};
}

/*! \relates NMatrix

  \brief Calculate a rotation matrix from a vector which encodes a
  rotation axis and angle.

  This is a right-handed expression, so all rotation axis are
  expected to be right handed.
*/
template <class T> inline NMatrix<T, 3> Rodrigues(const NVector<T, 3> &v) {
  const T theta = v.nrm();
  if (theta == 0)
    return NMatrix<T, 3>::identity();

  const NVector<T, 3> axis = v / theta;

  return NMatrix<T, 3>::identity() + std::sin(theta) * cross(axis) +
         (1 - std::cos(theta)) *
             (Dyadic(axis, axis) - NMatrix<T, 3>::identity());
}

/*! \relates NMatrix
  \brief Calculate the eigenvectors and values of a 3x3 NMatrix.
*/
inline std::pair<std::array<NVector<double, 3>, 3>, std::array<double, 3>>
symmetric_eigen_decomposition(const NMatrix<double, 3> &M) {
#ifdef MAGNET_DEBUG
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (M(i, j) != M(j, i))
        M_throw() << "Cannot perform an eigen decomposition of a matrix which "
                     "is not symmetric using this function!";
#endif

  double V[3][3];
  double d[3];
  double e[3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      V[i][j] = M(i, j);

  detail::tred2(V, d, e);
  detail::tql2(V, d, e);

  std::array<double, 3> eigenvals = {{d[0], d[1], d[2]}};
  std::array<NVector<double, 3>, 3> eigenvecs = {
      {NVector<double, 3>{V[0][0], V[1][0], V[2][0]},
       NVector<double, 3>{V[0][1], V[1][1], V[2][1]},
       NVector<double, 3>{V[0][2], V[1][2], V[2][2]}}};
  // Now output the eigenvectors and eigenvalues
  return std::make_pair(eigenvecs, eigenvals);
}

/*! \relates NVector
  \name NMatrix input/output operators
  \{
*/

/*! \brief XML output operator. */
template <class T, size_t N>
inline magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                          const NMatrix<T, N> &A) {
  char name[2] = "x";

  for (size_t iDim = 0; iDim < N; ++iDim) {
    name[0] = 'x' + iDim;
    XML << magnet::xml::tag(name);

    for (size_t jDim = 0; jDim < N; ++jDim) {
      char name2[2] = "x";
      name2[0] = 'x' + jDim;
      XML << magnet::xml::attr(name2) << A(iDim, jDim);
    }

    XML << magnet::xml::endtag(name);
  }

  return XML;
}

/*! \brief XML input operator. */
template <class T, size_t N>
NMatrix<T, N> &operator<<(NMatrix<T, N> &data, const magnet::xml::Node &XML) {
  char name[2] = "x";
  for (size_t iDim = 0; iDim < NDIM; ++iDim) {
    name[0] = 'x' + iDim;

    for (size_t jDim = 0; jDim < NDIM; ++jDim) {
      char name2[2] = "x";
      name2[0] = 'x' + jDim;

      data(iDim, jDim) = XML.getNode(name).getAttribute(name2).as<double>();
    }
  }
  return data;
}
/*! \} */

/*! \relates NVector
  \name Elementwise operations on NVectors
  \{
 */

template <class T1, class T2, size_t N>
inline NMatrix<decltype(T1() * T2()), N>
elementwiseMultiply(const NMatrix<T1, N> &A, const NMatrix<T2, N> &B) {
  NMatrix<decltype(T1() * T2()), N> retval;
  for (size_t i(0); i < 3; ++i)
    for (size_t j(0); j < 3; ++j)
      retval(i, j) = A(i, j) * B(i, j);
  return retval;
}

template <class T, size_t N>
inline NMatrix<T, N> elementwiseMin(const NMatrix<T, N> &A,
                                    const NMatrix<T, N> &B) {
  NMatrix<T, N> retval;
  for (size_t i(0); i < 3; ++i)
    for (size_t j(0); j < 3; ++j)
      retval(i, j) = std::min(A(i, j), B(i, j));
  return retval;
}

template <class T, size_t N>
inline NMatrix<T, N> elementwiseMax(const NMatrix<T, N> &A,
                                    const NMatrix<T, N> &B) {
  NMatrix<T, N> retval;
  for (size_t i(0); i < 3; ++i)
    for (size_t j(0); j < 3; ++j)
      retval(i, j) = std::max(A(i, j), B(i, j));
  return retval;
}

/*! \} */

typedef NMatrix<> Matrix;
} // namespace math
} // namespace magnet

namespace coil {
typedef ::magnet::math::NMatrix<> Matrix;
}
namespace dynamo {
typedef ::magnet::math::NMatrix<> Matrix;
}

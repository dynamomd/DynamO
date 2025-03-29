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
#include <magnet/math/bisect.hpp>
#include <magnet/math/symbolic.hpp>

namespace magnet {
namespace intersection {
/*! \brief Implementation of \ref next_root which utilises \ref
  solve_real_roots or \ref solve_real_positive_roots_poly.
*/
template <class Real, size_t Order, char Letter>
Real next_negative(const ::magnet::math::Polynomial<Order, Real, Letter> &f,
                   double t_origin = 0) {
  using namespace magnet::math;
  double first_root = HUGE_VAL;
  double second_root = HUGE_VAL;
  containers::StackVector<Real, Order> roots;
  typename containers::StackVector<Real, Order>::const_iterator it;

  // Where the equation is solved by radicals (3rd or lower order
  // polynomials), use those solutions. These return negative AND
  // positive roots so filter the negative roots out.
  if (Order < 4) {
    roots = solve_real_roots(f);
    it = roots.begin();
    // skip negative roots
    while ((it != roots.end()) && (*it < t_origin))
      ++it;
  } else {
    // For higher order polynomials, determine the positive roots
    roots = solve_real_positive_roots_poly<
        PolyRootBounder::VAS, PolyRootBisector::TOMS748, Order, Real, Letter>(
        f);
    it = roots.begin();
  }

  // Grab the next two positive roots (if they exist)
  if (it != roots.end()) {
    first_root = *it;
    if ((++it) != roots.end())
      second_root = *it;
  }

  // There are two numerical cases to be avoided.
  //
  // First, although f(0)>0 on entry to this function, there may be
  // a sign change numerically at x=0 due to finite precision. The
  // only way to test for this is to actually test if an unexpected
  // sign change has happened sometime before the next root (or
  // infinity if no root has happened). If there is no root change,
  // then x=1 is a logical point to sample the sign of the
  // root. Otherwise we must sample in-between now and the next
  // root.
  const double sample_point =
      (first_root == HUGE_VAL) ? t_origin + 1.0 : (t_origin + first_root) / 2;
  if (eval(f, Variable<Letter>() == sample_point) < 0)
    // They're already approaching, return an instant collision
    return t_origin;

  // Second, the root detected in the future may have an odd or
  // even multiplicity. Odd roots are actual sign transitions
  // whereas even roots are not. Therefore it is not enough to
  // simply detect the root, we must test that it actually changes
  // sign. Numerically determining the multiplicity of a root of a
  // floating point polynomial (or other function) is futile. We
  // only care if a sign change occurs.

  if (first_root != HUGE_VAL)
    do {
      const double sample_point = (second_root == HUGE_VAL)
                                      ? first_root + 1.0
                                      : (first_root + second_root) / 2;
      if (eval(f, Variable<Letter>() == sample_point) < 0)
        return first_root;
      first_root = second_root;
      second_root = ((++it) != roots.end()) ? *it : HUGE_VAL;
    } while (first_root != HUGE_VAL);

  return HUGE_VAL;
}

/*!  \brief Implementation of the generic stable event-detection
  algorithm.

  For this generic implementation to work, the function f must
  have implementations of the functions \ref
  magnet::math::derivative, \ref magnet::math::shift_function, and
  \ref magnet::math::next_root defined.
*/
template <class Function> double nextEvent(Function f, double t_origin = 0) {
  ::magnet::math::Variable<'t'> t;

  // Determine the derivative
  const auto df = ::magnet::math::derivative(f, t);

  const double f_start = eval(f, t == t_origin);
  const double df_start = eval(df, t == t_origin);

  // Check if starting was not overlapped.
  if (f_start >= 0) {
    return next_negative(f, t_origin);
  }

  // If we're approaching then the current time is the time of
  // the next event
  if (df_start < 0)
    return 0;

  // We need to find when the derivative next turns negative.
  const double df_next_root = next_negative(df, t_origin);

  // If the overlap function never becomes negative, then there is
  // never an event.
  if (df_next_root == HUGE_VAL)
    return HUGE_VAL;

  // If it turns around while still overlapped/in contact, then
  // the turning point is the next event.
  if (eval(f, t == df_next_root) <= 0)
    return df_next_root;

  // Ok, try searching after the maxima for events
  return next_negative(f, t_origin + df_next_root);
}

/*! \brief Calculate the interval until the 1st order Polynomial
    is negative and has a negative derivative.

    This is a specialisation for linear polynomials to close the
    recursive definition of the general nextEvent implementation.

    \param f The Polynomial under consideration.
*/
template <char Var>
inline double nextEvent(const ::magnet::math::Polynomial<1, double, Var> &f) {
  // If the gradient is not negative now, it never will be
  if (f[1] >= 0)
    return HUGE_VAL;
  // Return the time of the root, or now if we're past the root
  return std::max(0.0, -f[0] / f[1]);
}

/*! \brief Calculate the interval until the 2nd order Polynomial
    is negative and has a negative derivative.

    This is an optimised case, as many interactions use quadratic
    overlap functions.

    \param f The Polynomial under consideration.
*/
template <char Var>
inline double nextEvent(const ::magnet::math::Polynomial<2, double, Var> &f) {
  // If the polynomial is linear, drop to that solution
  if (f[2] == 0)
    return nextEvent(::magnet::math::change_order<1>(f));

  const double arg = f[1] * f[1] - 4 * f[2] * f[0];

  if (f[2] < 0) {
    // Polynomial limits towards overlap at t -> +inf

    // If there are no roots, it never escapes overlap, return the
    // time of the turning point or now if it is in the past
    if (arg <= 0)
      return std::max(0.0, -f[1] / (2 * f[2]));

    // There are roots. select a stable form of the quadratic to
    // compute the largest root.
    if (f[1] > 0)
      return std::max(0.0, (-f[1] - std::sqrt(arg)) / (2 * f[2]));
    else
      return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
  }
  /* else (f[2] > 0)*/

  // Interactions only happen if there are roots and we're in the
  // region between the first root and the turning point
  if ((f[1] >= 0) || arg <= 0)
    return HUGE_VAL;

  // Return the time of the root using a stable quadratic formula
  return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
}
} // namespace intersection
} // namespace magnet

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
#include <array>
#include <limits>
#include <magnet/math/bisect.hpp>
#include <magnet/math/frenkelroot.hpp>

namespace magnet {
namespace intersection {
namespace detail {
/*! \brief A class which takes the derivative of an overlap
  function.

  The derivative is taken by reference, so the original overlap
  function must not fall out of scope.
 */
template <class Base, int derivative = 1> class FDerivative {
  const Base &_f;

public:
  FDerivative(const Base &f) : _f(f) {}

  template <size_t first_deriv = 0, size_t nderivs = 1>
  std::array<double, nderivs> eval(double dt) const {
    return _f.template eval<first_deriv + derivative, nderivs>(dt);
  }
  template <size_t d> double max() const {
    return _f.template max<d + derivative>();
  }
};

template <class Base> class FBisect_Wrapper {
  const Base &_f;

public:
  FBisect_Wrapper(const Base &f) : _f(f) {}
  double operator()(const double dt) const {
    return _f.template eval<0>(dt).front();
  }
};

/*! \brief A numerical root finder based on Halley's method.

  This is losely based around the boost implementation of
  Halley's method which includes additional strategies to
  attempt to ensure convergence. We cannot use the boost method
  directly as it seems to assume a root is definitely contained
  within the bounds whereas we are speculatively searching for a
  root.

  This implementation attempts a Halley's method step by
  default, but if this has numerical difficulties, it switches
  to a Newton-Raphson step.

  The search window [t_min,t_max] is updated with each step to
  make sure that the steps are converging (delta is reducing
  with each step).
  If this does not happen, a step will head out of bounds
  possibly indicating a) that our initial guess is too far from
  a root for the quadratic approximation to work or b) a
  boundary is very near the root. To deal with case b), the
  mid-point between the current location and the boundary in the
  direction of the step (t_min/t_max) is used as the next
  step. This is a last-ditch attempt to get the method to
  converge.

  The overall key to success is an appropriate selection of the
  binary_digits argument. If we try to solve roots to a higher
  precision than f allows, we might enter endless cycles due to
  precision, therefore binary_digits cannot be too large.
 */
template <class F>
std::pair<bool, double> halleySearch(const F &f, double t_guess, double t_min,
                                     double t_max, const double fprecision,
                                     const int binary_digits,
                                     size_t iterations) {
  const double digitfactor = std::ldexp(1.0, 1 - binary_digits);

  // This is recalculated at the end of each step to allow a check for
  // convergence.
  auto fval = f.template eval<0, 3>(t_guess);
  do {
    // Check if we're at a root already
    if (fval[0] == 0)
      return std::pair<bool, double>(true, t_guess);

    // If we have zero derivatives, we can't make progress so
    // abort to allow the algorithm to try again from somewhere
    // else
    if ((fval[1] == 0) && (fval[2] == 0))
      break;

    // Calculate the numerator and denominator terms of Halley's method.
    const double denom = 2 * fval[0];
    const double num = 2 * fval[1] - fval[0] * (fval[2] / fval[1]);

    // Calculate the delta, taking care over the evaluation and
    // switching to a Newton step if required.
    double delta;
    if ( // Check that we had a second derivative
        (fval[2] == 0)
        // Check for overflow
        || ((std::abs(num) < 1) &&
            (std::abs(denom) >=
             std::abs(num) * std::numeric_limits<double>::max())))
      delta = -fval[0] / fval[1];
    else {
      // Perform a Halley step
      delta = -denom / num;
      // Check that both Newton and Halley methods agree on the
      // direction, if they don't it could be a cancellation
      // error. We then switch to a Newton step, with a fixed
      // maximum magnitude
      if (-delta * fval[1] / fval[0] < 0) {
        delta = -fval[0] / fval[1];
        const double max_delta = 2 * std::abs(t_guess);
        if (std::abs(delta) > max_delta)
          delta = std::copysign(max_delta, delta);
      }
    }

    // Calculate the new step
    double t_new_guess = t_guess + delta;

    // Check we've not gone out of range
    if ((t_new_guess < t_min) || (t_new_guess > t_max))
      break;

    ////Accept the step and update the bounds to the old guesses.
    if (t_new_guess > t_guess)
      t_min = t_guess;
    else
      t_max = t_guess;
    // Then update the current guess.
    t_guess = t_new_guess;

    // Check if we've converged (and calculate the next iteration's data)
    fval = f.template eval<0, 3>(t_guess);
    if (std::abs(t_guess * digitfactor) >= std::abs(delta) &&
        (std::abs(fval[0]) < fprecision))
      return std::pair<bool, double>(true, t_guess);
  } while (--iterations);

  // Failed! We've not found a valid root.
  return std::pair<bool, double>(false, HUGE_VAL);
}

template <class F>
std::pair<bool, double>
nextDecreasingRoot(const F &f, double t_min, double t_max,
                   double fprecision = 1e-10,
                   size_t restarts = std::numeric_limits<size_t>::max() - 1,
                   const size_t halley_binary_digits = 32,
                   const size_t halley_iterations = 500) {
  // Make things clearer using enums for high/low boundary marking
  enum { LOW = 0, HIGH = 1 };
  int active_boundary = LOW;

  // Cache this value as it may be expensive to compute
  const double f2max = f.template max<2>();

  // Loop while we still have a valid search window and
  // are not limiting the number of restarts

  // Store the initial sign of the function at t_min and
  // t_max. This is used to combat precision errors which cause
  // the boundaries to pass over a root when updated.
  const bool t_min_sign = std::signbit(f.template eval<0>(t_min).front());
  const bool t_max_sign = std::signbit(f.template eval<0>(t_max).front());
  double old_t_min = t_min, old_t_max = t_max;

  ++restarts;
  while ((t_min < t_max) && (--restarts)) {
    double &t_current = (active_boundary == HIGH) ? t_max : t_min;
    const auto fval = f.template eval<0, 2>(t_current);

    // Improve the boundary. The copysign should guarantee that
    // we have roots either side of t_current.
    auto boundary_roots = math::quadraticEquation(
        -0.5 * std::copysign(f2max, fval[0]), fval[1], fval[0]);

    // Check that the last update didn't cause a sign change on
    // at the boundary (indicating it skipped over a root). If
    // it did, we can bisect a root here!
    if (active_boundary == HIGH) {
      if (std::signbit(fval[0]) == t_max_sign) {
        old_t_max = t_max;
        t_max += std::min(boundary_roots.first, boundary_roots.second);
      } else
        return std::pair<bool, double>(
            true, magnet::math::bisect(FBisect_Wrapper<F>(f), old_t_max, t_max,
                                       fprecision));
    } else {
      if (std::signbit(fval[0]) == t_min_sign) {
        old_t_min = t_min;
        t_min += std::max(boundary_roots.first, boundary_roots.second);
      } else
        return std::pair<bool, double>(
            true, magnet::math::bisect(FBisect_Wrapper<F>(f), old_t_min, t_min,
                                       fprecision));
    }

    // Switch the working boundary for the next iteration (if
    // we have a finite upper bound, otherwise stick with the
    // lower bound)
    if (t_max != HUGE_VAL)
      active_boundary = !active_boundary;

    // Now search for a root using halley's method, starting
    // from the current boundary
    auto search_result = halleySearch(f, t_current, t_min, t_max, fprecision,
                                      halley_binary_digits, halley_iterations);

    // If searching failed, restart from the other bound (and
    // update it)
    if (!search_result.first)
      continue;

    {
      // Searching succeeded check for earlier roots in this
      // interval using recursion.
      const double current_root = search_result.second;
      const double f1 = f.template eval<1>(current_root).front();
      const double inner_t_max = current_root - 2.0 * std::abs(f1 / f2max);
      auto check_result = nextDecreasingRoot(
          f, t_min, inner_t_max, halley_binary_digits, halley_iterations);
      if (check_result.second != HUGE_VAL)
        // We have found an earlier root. Return this one instead.
        return check_result;

      // There are no earlier roots, current_root is the
      // next root in the region. Return it if its
      // approaching
      if (f1 < 0)
        return std::pair<bool, double>(true, current_root);

      // Use it as a new lower bound if its a receeding
      // root and carry on.
      t_min = current_root + 2.0 * std::abs(f1 / f2max);
    }
  }

  // We've failed to find a root
  if (restarts == 0)
    return std::pair<bool, double>(false, t_min);
  else
    return std::pair<bool, double>(true, HUGE_VAL);
}
} // namespace detail

/*! \brief A generic implementation of the stable EDMD algorithm.
  \tparam T The type of the overlap function which is being solved.
  \param f The overlap function.
  \param t_min The minimum time to start the search from.
  \param t_max The maximum time to start the search from (may be HUGE_VAL).
*/
template <class T>
std::pair<bool, double> nextEvent(const T &f, const double t_min,
                                  const double t_max) {
  const auto fval = f.template eval<0, 2>(t_min);

  // If particles are not in contact, just search for the next contact
  if (fval[0] > 0)
    return detail::nextDecreasingRoot(f, t_min, t_max);

  // Particles are either in contact or overlapped. Check if they're approaching
  if (fval[1] < 0)
    return std::pair<bool, double>(true, 0.0);

  // Overlapping but they're moving away from each other. Determine
  // when they reach their next maximum separation (it may be the
  // current time if f1==0).
  detail::FDerivative<T> fprime(f);
  std::pair<bool, double> derivroot =
      detail::nextDecreasingRoot(fprime, t_min, t_max);
  // Check if they just keep retreating from each other, which means that they
  // never interact
  if (derivroot.second == HUGE_VAL)
    return std::pair<bool, double>(false, HUGE_VAL);

  // If they are still overlapping at this time, it doesn't matter
  // if derivroot is a virtual (recalculate) event or an actual
  // turning point. We can just return it and either have a
  // collision then or recalculate then.
  if (f.template eval<0>(derivroot.second).front() < 0)
    return derivroot;

  // Real or virtual, the derivroot contains a time before the
  // next interaction which is outside the invalid state, we just
  // use this as our lower bound and carry on the search.
  return detail::nextDecreasingRoot(f, derivroot.second, t_max);
}
} // namespace intersection
} // namespace magnet

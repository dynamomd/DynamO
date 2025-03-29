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
#include <limits>
#include <magnet/math/quadratic.hpp>

/*
   This work is heavily derived from the public domain work of Don
   Herbison-Evans. The original code is available in
   src/magnet/test/quartic_original.cpp. The code has been refactored
   to change its coding style. Any changes to the function are listed
   below.

   Oct 2013: Following an email from Florian Bruckner, The special
   case of ((p == 0) && (q == 0)) has been modified to return 1 or 0
   roots instead of 3, and the std::cbrt function is used instead of
   pow.
*/

namespace magnet {
namespace math {
namespace detail {

/*! \brief Uses a quadratic scheme to polish up a root,
    switching to a bisection scheme if it manages to bracket the
    root.
 */
inline void cubicNewtonRootPolish(const double &p, const double &q,
                                  const double &r, double &root) {
  // Stored for comparison later
  double error = ((root + p) * root + q) * root + r;
  const size_t maxiterations = 2;
  for (size_t it = 0; (it < maxiterations) && (error != 0); ++it) {
    double deriv = (3 * root + 2 * p) * root + q;
    double dderiv = 6 * root + 2 * p;

    // Try to use a quadratic scheme to improve the root
    std::pair<double, double> roots;
    try {
      roots = quadraticEquation(0.5 * dderiv, deriv, error);
      root += (std::abs(roots.first) < std::abs(roots.second)) ? roots.first
                                                               : roots.second;
    } catch (NoQuadraticRoots &) {
      // Switch to a linear scheme if the quadratic fails
      if (deriv == 0)
        return;
      root -= error / deriv;
    }

    error = ((root + p) * root + q) * root + r;
  }
}
} // namespace detail

// Please read  http://linus.it.uts.edu.au/~don/pubs/solving.html
// For solving cubics like x^3 + p * x^2 + q * x + r == 0
inline size_t cubicSolve(const double &p, const double &q, const double &r,
                         double &root1, double &root2, double &root3) {
  static const double maxSqrt = std::sqrt(std::numeric_limits<double>::max());

  if (r == 0) {
    // no constant term, so divide by x and the result is a
    // quadratic, but we must include the trivial x = 0 root
    try {
      std::pair<double, double> roots = quadraticEquation(1.0, p, q);
      root1 = roots.first;
      root2 = roots.second;
      root3 = 0;
      if (root1 < root2)
        std::swap(root1, root2);
      if (root2 < 0) {
        std::swap(root2, root3);
        if (root1 < 0)
          std::swap(root1, root2);
      }
      return 3;
    } catch (NoQuadraticRoots &) {
    }

    root1 = 0;
    return 1;
  }

  if ((p == 0) && (q == 0)) {
    // Special case
    // Equation is x^3 == -r
    if (r > 0)
      return 0;
    root1 = std::cbrt(-r);
    return 1;
  }

  if ((p > maxSqrt) || (p < -maxSqrt)) {
    // Equation limits to x^3 + p * x^2 == 0
    root1 = -p;
    return 1;
  }

  if (q > maxSqrt) {
    // Special case, if q is large and the root is -r/q,
    // The x^3 term is negligble, and all other terms cancel.
    root1 = -r / q;
    return 1;
  }

  if (q < -maxSqrt) {
    // Special case, equation is x^3 + q x == 0
    root1 = -std::sqrt(-q);
    return 1;
  }

  if ((r > maxSqrt) || (r < -maxSqrt)) {
    // Another special case
    // Equation is x^3 == -r
    root1 = -std::cbrt(r);
    return 1;
  }

  double v = r + (2.0 * p * p / 9.0 - q) * (p / 3.0);

  if ((v > maxSqrt) || (v < -maxSqrt)) {
    root1 = -p;
    return 1;
  }

  double uo3 = q / 3.0 - p * p / 9.0;
  double u2o3 = uo3 + uo3;

  if ((u2o3 > maxSqrt) || (u2o3 < -maxSqrt)) {
    if (p == 0) {
      if (q > 0) {
        root1 = -r / q;
        return 1;
      }

      if (q < 0) {
        root1 = -std::sqrt(-q);
        return 1;
      }

      root1 = 0;
      return 1;
    }

    root1 = -q / p;
    return 1;
  }

  double uo3sq4 = u2o3 * u2o3;
  if (uo3sq4 > maxSqrt) {
    if (p == 0) {
      if (q > 0) {
        root1 = -r / q;
        return 1;
      }

      if (q < 0) {
        root1 = -std::sqrt(-q);
        return 1;
      }

      root1 = 0;
      return 1;
    }

    root1 = -q / p;
    return 1;
  }

  double j = (uo3sq4 * uo3) + v * v;

  if (j > 0) { // Only one root (but this test can be wrong due to a
    // catastrophic cancellation in j
    //(i.e. (uo3sq4 * uo3) == v * v)
    double w = std::sqrt(j);
    if (v < 0)
      root1 =
          std::cbrt(0.5 * (w - v)) - (uo3)*std::cbrt(2.0 / (w - v)) - p / 3.0;
    else
      root1 =
          uo3 * std::cbrt(2.0 / (w + v)) - std::cbrt(0.5 * (w + v)) - p / 3.0;

    // We now polish the root up before we use it in other calculations
    detail::cubicNewtonRootPolish(p, q, r, root1);

    // We double check that there are no more roots by using a
    // quadratic formula on the factored problem, this helps when
    // the j test is wrong due to numerical error.

    // We have a choice of either -r/root1, or q -
    //(p+root1)*root1 for the constant term of the quadratic.
    //
    // The division one usually results in more accurate roots
    // when it finds them but fails to detect real roots more
    // often than the multiply.
    try {
      std::pair<double, double> roots =
          quadraticEquation(1.0, p + root1, -r / root1);
      root2 = roots.first;
      root3 = roots.second;
      return 3;
    } catch (NoQuadraticRoots &) {
    }

    // However, the multiply detects roots where there are none,
    // the division does not. So we must either accept missing
    // roots or extra roots, here we choose missing roots
    //
    // if (quadSolve(q-(p+root1)*root1, p + root1, 1.0, root2, root3))
    //   return 3;

    return 1;
  }

  if (uo3 >= 0) { // Multiple root detected
    root1 = root2 = root3 = std::cbrt(v) - p / 3.0;
    return 3;
  }

  double muo3 = -uo3;
  double s;
  if (muo3 > 0) {
    s = std::sqrt(muo3);
    if (p > 0)
      s = -s;
  } else
    s = 0;

  double scube = s * muo3;
  if (scube == 0) {
    root1 = -p / 3.0;
    return 1;
  }

  double t = -v / (scube + scube);
  double k = std::acos(t) / 3.0;
  double cosk = std::cos(k);
  root1 = (s + s) * cosk - p / 3.0;

  double sinsqk = 1.0 - cosk * cosk;
  if (sinsqk < 0)
    return 1;

  double rt3sink = std::sqrt(3.0) * std::sqrt(sinsqk);
  root2 = s * (-cosk + rt3sink) - p / 3.0;
  root3 = s * (-cosk - rt3sink) - p / 3.0;

  detail::cubicNewtonRootPolish(p, q, r, root1);
  detail::cubicNewtonRootPolish(p, q, r, root2);
  detail::cubicNewtonRootPolish(p, q, r, root3);

  return 3;
}
} // namespace math
} // namespace magnet

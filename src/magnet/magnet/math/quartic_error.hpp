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

#include <algorithm>
#include <array>
#include <cmath>

/*
   This work is heavily derived from the public domain work of Don
   Herbison-Evans. The original code is available in
   src/magnet/test/quartic_original.cpp. The code has been refactored
   to change its coding style. Any changes to the function are listed
   below.
*/

namespace magnet {
namespace math {
inline double quarticError(const double &a, const double &b, const double &c,
                           const double &d, const double roots[4],
                           const size_t rootCount) {
  std::array<double, 4> errors = {{0, 0, 0, 0}};

  for (size_t root = 0; root < rootCount; ++root) {
    const double value =
        (((roots[root] + a) * roots[root] + b) * roots[root] + c) *
            roots[root] +
        d;

    if (value == 0) {
      errors[root] = 0;
      continue;
    }

    const double deriv =
        ((4 * roots[root] + 3 * a) * roots[root] + 2 * b) * roots[root] + c;

    if (deriv != 0)
      errors[root] = std::abs(value / deriv);
    else {
      const double secDeriv = (12 * roots[root] + 6 * a) * roots[root] + 2 * b;
      if (secDeriv != 0)
        errors[root] = std::sqrt(std::abs(value / secDeriv));
      else {
        const double thirdDeriv = 24 * roots[root] + 6 * a;
        if (thirdDeriv != 0)
          errors[root] = std::pow(std::abs(value / thirdDeriv), 1.0 / 3.0);
        else
          errors[root] = std::sqrt(std::sqrt(std::abs(value) / 24));
      }
    }
  }

  return *std::max_element(errors.begin(), errors.begin() + rootCount);
}
} // namespace math
} // namespace magnet

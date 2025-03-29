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
#include <boost/circular_buffer.hpp>
#include <exception>
#include <magnet/math/vector.hpp>
#include <tuple>
#include <utility>
#include <vector>

namespace magnet {
namespace math {
/*! \brief This class facilitates the exact averaging of
    properties during an event driven simulation.

    As many simulation properties are constant between events, we
    can actually calculate exact statistics on the fluctuations of
    these properties. This class calculates the mean, square mean,
    minimum and maximum of a given property.
 */
template <class T> class TimeAveragedProperty {
public:
  TimeAveragedProperty()
      : _current_value(T()), _zero_moment(0), _first_moment(T()),
        _second_moment(T()), _min(T()), _max(T()) {}

  void swapAverages(TimeAveragedProperty &op) {
    std::swap(_zero_moment, op._zero_moment);
    std::swap(_first_moment, op._first_moment);
    std::swap(_second_moment, op._second_moment);
    std::swap(_min, op._min);
    std::swap(_max, op._max);
  }

  void init(const T &value) { _current_value = _min = _max = value; }

  TimeAveragedProperty &operator=(const T &value) {
    _current_value = value;
    _min = elementwiseMin(_current_value, _min);
    _max = elementwiseMax(_current_value, _max);
    return *this;
  }

  TimeAveragedProperty &operator+=(const T &change) {
    return operator=(_current_value + change);
  }

  void stream(double dt) {
    _zero_moment += dt;
    _first_moment += dt * _current_value;
    _second_moment += dt * elementwiseMultiply(_current_value, _current_value);
  }

  double time() const { return _zero_moment; }

  T mean() const {
    return (_zero_moment == 0) ? current() : (_first_moment / _zero_moment);
  }

  T meanSqr() const {
    return (_zero_moment == 0) ? (current() * current())
                               : (_second_moment / _zero_moment);
  }

  T min() const { return _min; }

  T max() const { return _max; }

  T current() const { return _current_value; }

protected:
  T _current_value;
  double _zero_moment;
  T _first_moment;
  T _second_moment;
  T _min;
  T _max;
};
} // namespace math
} // namespace magnet

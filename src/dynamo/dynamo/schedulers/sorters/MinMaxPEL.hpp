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
#include <dynamo/eventtypes.hpp>
#include <magnet/containers/MinMaxHeap.hpp>
#include <string>

namespace dynamo {
/*! A MinMax heap used for Particle Event Lists

  There is a trick used here to speed up comparisons between
  MinMaxHeaps.  The top element is set to
  std::numeric_limits<float>::infinity(), whenever the queue is cleared, or
  pop'd empty. This means no conditional logic is required to deal with the
  comparison of empty queues.
*/
template <size_t Size> class MinMaxPEL {
  magnet::containers::MinMaxHeap<Event, Size> _store;

public:
  static const bool partial_invalidate_support = false;

  MinMaxPEL() { clear(); }

  inline void push(const Event &e) {
    if (!_store.full())
      _store.insert(e);
    else {
      if (e < _store.bottom())
        _store.replaceMax(e);
      _store.unsafe_bottom()._type = RECALCULATE;
      _store.unsafe_bottom()._source = SCHEDULER;
    }
  }

  inline void clear() {
    _store.clear();
    (*_store.begin()) = Event();
  }

  inline size_t size() const { return _store.size(); }

  inline bool empty() const { return _store.empty(); }

  inline void pop() {
    _store.pop();
    if (_store.empty())
      clear();
  }

  inline Event top() const { return *_store.begin(); }

  inline bool operator>(const MinMaxPEL &o) const { return top() > o.top(); }

  inline bool operator<(const MinMaxPEL &o) const { return top() < o.top(); }

  inline void stream(const double dt) {
    for (Event &event : _store)
      event._dt -= dt;
  }

  inline void rescaleTimes(const double scale) {
    for (Event &event : _store)
      event._dt *= scale;
  }

  inline void swap(MinMaxPEL &rhs) { _store.swap(rhs._store); }

  static inline std::string name() { return "MinMax" + std::to_string(Size); }
};
} // namespace dynamo

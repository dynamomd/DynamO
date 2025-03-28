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
#include <dynamo/eventtypes.hpp>
#include <functional>
#include <vector>

namespace dynamo {
class HeapPEL {
  std::vector<Event> _store;

public:
  static const bool partial_invalidate_support = false;

  inline void push(Event e) {
    _store.push_back(e);
    std::push_heap(_store.begin(), _store.end(), std::greater<Event>());
  }

  inline void clear() { _store.clear(); }

  inline size_t size() const { return _store.size(); }

  inline bool empty() const { return _store.empty(); }

  inline void pop() {
    std::pop_heap(_store.begin(), _store.end(), std::greater<Event>());
    _store.pop_back();
  }

  inline Event top() const {
    if (!empty())
      return _store.front();
    else
      return Event();
  }

  inline bool operator>(const HeapPEL &FEL) const { return top() > FEL.top(); }

  inline bool operator<(const HeapPEL &FEL) const { return top() < FEL.top(); }

  inline void stream(const double dt) {
    for (Event &event : _store)
      event._dt -= dt;
  }

  inline void rescaleTimes(const double scale) {
    for (Event &event : _store)
      event._dt *= scale;
  }

  inline void swap(HeapPEL &rhs) { std::swap(_store, rhs._store); }

  static inline std::string name() { return "Heap"; }
};
} // namespace dynamo

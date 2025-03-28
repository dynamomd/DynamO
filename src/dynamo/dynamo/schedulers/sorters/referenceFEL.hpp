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
#include <dynamo/schedulers/sorters/FEL.hpp>

namespace dynamo {
/*! \brief A slow, but exact FEL implementation.

  This is used as a reference
 */
class ReferenceFEL : public FEL {
  std::vector<Event> _store;

public:
  virtual void init(const size_t N) {}

  virtual void clear() { _store.clear(); }
  virtual bool empty() { return _store.empty(); }

  virtual void push(Event e) { _store.push_back(e); }

  virtual void rescaleTimes(const double f) {
    for (Event &e : _store)
      e._dt *= f;
  }

  virtual void stream(const double dt) {
    for (Event &e : _store)
      e._dt -= dt;
  }

  virtual Event top() {
    if (empty())
      M_throw() << "Event queue is empty!";
    return *std::min_element(_store.begin(), _store.end());
  }

  virtual void invalidate(const size_t id) {
    auto test = [=](const Event &e) {
      return (e._particle1ID == id) ||
             ((e._source == INTERACTION) && (e._particle2ID == id));
    };

    _store.erase(std::remove_if(_store.begin(), _store.end(), test),
                 _store.end());
  }

  virtual void pop() {
    _store.erase(std::min_element(_store.begin(), _store.end()));
  }

private:
  virtual void outputXML(magnet::xml::XmlStream &) const {}
};
} // namespace dynamo

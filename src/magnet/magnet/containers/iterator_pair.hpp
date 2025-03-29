/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2008  Todd Wease <->

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Based on a an original implementation kindly released by Todd
    Wease.
  */
#pragma once

namespace magnet {
namespace containers {
/*! \brief This container allows a pair of iterators to be used in
    a range-based for loop.

    An example usage, using the helper function, is this:

    \code
    std::vector<int> vector;
    for (auto& value : IteratorPairRange<>)
    \endcode
*/
template <class Iterator> class IteratorPairRange {
public:
  IteratorPairRange(Iterator begin, Iterator end) : _begin(begin), _end(end) {}
  const Iterator &begin() const { return _begin; }
  const Iterator &end() const { return _end; }

private:
  const Iterator _begin;
  const Iterator _end;
};

/*! \brief Helper function for IteratorPairRange. */
template <class Iterator>
IteratorPairRange<Iterator> makeIteratorRange(Iterator begin, Iterator end) {
  return IteratorPairRange<Iterator>(begin, end);
}

} // namespace containers
} // namespace magnet

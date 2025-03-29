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
#include <iomanip>
#include <magnet/string/searchreplace.hpp>

namespace magnet {
namespace string {
namespace detail {
/*! \brief Class to track the number of the line.
 */
class LineNum {
public:
  LineNum(size_t count_width) : _count(1), _width(count_width) {}

  LineNum &operator++() {
    ++_count;
    return *this;
  }

  friend std::ostream &operator<<(std::ostream &os, LineNum &id) {
    return os << std::setw(id._width) << id._count++ << ": ";
  }

protected:
  size_t _count;
  size_t _width;
};
} // namespace detail
/*! \brief Formats text by adding line numbers.

  \param in The string containing the source text.

  \returns A formatted version of the text.
 */
inline std::string add_line_numbers(std::string in) {
  std::ostringstream os;
  size_t lines = std::count(in.begin(), in.end(), '\n') + 1;
  size_t digits = 0;
  while (lines) {
    lines /= 10;
    ++digits;
  }

  detail::LineNum number(digits);

  os << number;

  for (std::string::const_iterator cPtr = in.begin(); cPtr != in.end();
       ++cPtr) {
    os << *cPtr;
    switch (*cPtr) {
    case '\n':
      os << number;
      break;
    }
  }

  return os.str();
}
} // namespace string
} // namespace magnet

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
#include <string>

namespace magnet {
namespace string {
/*! \brief Search and replace elements in a std::string.
 * \param in The string to search within.
 * \param from A string giving the text sequence to replace
 * \param to A string with the replacement text sequence.
 * \returns The string "in" with all occurences of "from" replaced with "to".
 */
inline std::string search_replace(std::string in, const std::string &from,
                                  const std::string &to) {
  if (!in.empty()) {
    std::string::size_type toLen = to.length();
    std::string::size_type frLen = from.length();
    std::string::size_type loc = 0;

    while (std::string::npos != (loc = in.find(from, loc))) {
      in.replace(loc, frLen, to);
      loc += toLen;

      if (loc >= in.length())
        break;
    }
  }
  return in;
}
} // namespace string
} // namespace magnet

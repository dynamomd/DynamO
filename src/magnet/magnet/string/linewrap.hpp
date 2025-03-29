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
#include <locale>
#include <string>

namespace magnet {
namespace string {
namespace {
inline size_t validCharCount(const std::string in, size_t pos1, size_t pos2) {
  std::string data(in.begin() + pos1, in.begin() + pos2);
  // Search and remove console commands
  for (size_t i(0); i < data.size(); ++i)
    if (data[i] == '\033')
      for (size_t j(i + 1); j < data.size(); ++j)
        if (data[j] == 'm') {
          data.erase(data.begin() + i, data.begin() + j + 1);
          break;
          i = 0;
        }

  return data.size();
}
} // namespace

/*! \brief Line wraps a string in a neat way.

  This function will wrap lines before they overflow

  \param in The string to be wrapped.
  \param linelength The length of the line in characters.
  \tparam hyphenate_long_words If true, any word longer than a line will be
  broken using hyphens
  \returns The wrapped form of the in string.
 */
template <bool hyphenate_long_words>
inline std::string linewrap(std::string in, size_t linelength) {
  // The position of the newline or space character before the
  // current word
  std::string::size_type linestart = -1;

  // The position of the newline or space character before the
  // current line
  std::string::size_type wordstart = -1;

  // Loop over the string, character by character
  for (std::string::size_type curchar = 0; curchar <= in.size(); ++curchar) {
    // Are we at the end of a word?
    if ((in[curchar] == ' ') || (curchar == in.size())) {
      // Will this word take the line length over its maximum?
      if (validCharCount(in, linestart + 1, curchar) > linelength) {
        // Check if we can insert a line break before the current word.
        if (linestart != wordstart) {
          in[wordstart] = '\n';
          linestart = wordstart;
        }

        // Now check if we need to break the current word up (and hyphenate it)
        while (validCharCount(in, linestart + 1, curchar) > linelength) {
          if (hyphenate_long_words)
            in.insert(linestart + linelength, std::string("\n-"));
          else
            in.insert(linestart + linelength, std::string("\n"));

          linestart = linestart + linelength;
          curchar += 1 + hyphenate_long_words;
        }
      }

      // Mark the start of the next word
      wordstart = curchar;
    } else if (in[curchar] == '\n')
      wordstart = linestart = curchar;
  }
  return in;
}
} // namespace string
} // namespace magnet

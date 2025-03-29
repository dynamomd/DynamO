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

namespace magnet {
namespace gtk {
/*! \brief This function is designed to be used in a Gtk::Entry's
 signal_changed() callback.  It enforces that the text is
 formatted in the style of a numeric value.

 It allows the first character to be a minus or plus sign, there
 to be a sequence of digits, possibly split by one decimal point.

 An example use is \code
 entry.signal_changed().connect(sigc::bind<Gtk::Entry*>(&magnet::Gtk::forceNumericEntry,
 &entry)); \endcode
*/
inline void forceNumericEntry(::Gtk::Entry *textfield) {
  std::string value = textfield->get_text();

  bool hasPoint = false;
  bool hasExponent = false;
  std::string::iterator iPtr = value.begin();
  if ((*iPtr == '-') || (*iPtr == '+'))
    ++iPtr;

  while (iPtr != value.end())
    if (std::isdigit(*iPtr))
      ++iPtr;
    else if ((*iPtr == '.') && (!hasPoint)) {
      ++iPtr;
      hasPoint = true;
    } else if ((*iPtr == 'e') && (!hasExponent) && (iPtr != value.begin()) &&
               std::isdigit(*(iPtr - 1))) {
      // If the last value was a digit we can have a single e
      // argument, but don't allow decimal exponents
      hasExponent = true;
      hasPoint = true;
      ++iPtr;
      // Eat the sign of the exponent
      if ((*iPtr == '-') || (*iPtr == '+'))
        ++iPtr;
    } else
      iPtr = value.erase(iPtr);

  if (value[0] == '.')
    value.erase(value.begin());

  textfield->set_text(value);
}
} // namespace gtk
} // namespace magnet

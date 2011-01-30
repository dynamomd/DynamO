/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
  namespace Gtk {
    //!This function is designed to be used in a Gtk::Entry's signal_changed() callback.
    // It enforces that the text is formatted in the style of a numeric value
    // E.g.  entry.signal_changed().connect(sigc::bind<Gtk::Entry&>(&magnet::Gtk::forceNumericEntry, entry));
    inline void forceNumericEntry(::Gtk::Entry& textfield)
    {
      std::string value = textfield.get_text();
      
      bool hasPoint = false;
      for (std::string::iterator iPtr = value.begin(); iPtr != value.end();)
	if (std::isdigit(*iPtr))
	  ++iPtr;
	else
	  if ((*iPtr == '.') && (!hasPoint))
	    { ++iPtr; hasPoint = true; }
	  else
	    iPtr = value.erase(iPtr);
      
      if (value[0] == '.') value.erase(value.begin());

      textfield.set_text(value);
    }
  }
}

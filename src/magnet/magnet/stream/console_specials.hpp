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

namespace magnet
{
  namespace console {
    //! \brief All text after this command are bold.
    inline std::string bold  () { return "\033[1m";}
    //! \brief All text after this command are underlined.
    inline std::string underline  () { return "\033[4m";}
    //! \brief All text after this command are struck through.
    inline std::string strikethrough  () { return "\033[9m";}
    
    //! \brief Set the text color to black.
    inline std::string black_fg  () { return "\033[0;30m";}
    //! \brief Set the text color to red.
    inline std::string red_fg    () { return "\033[0;31m"; }
    //! \brief Set the text color to green.
    inline std::string green_fg  () { return "\033[0;32m"; }
    //! \brief Set the text color to yellow.
    inline std::string yellow_fg () { return "\033[0;33m"; }
    //! \brief Set the text color to blue.
    inline std::string blue_fg   () { return "\033[0;34m"; }
    //! \brief Set the text color to purple.
    inline std::string purple_fg () { return "\033[0;35m"; }
    //! \brief Set the text color to cyan.
    inline std::string cyan_fg   () { return "\033[0;36m"; }
    //! \brief Set the text color to white.
    inline std::string white_fg  () { return "\033[0;37m"; }

    //! \brief Set the background color to black.
    inline std::string black_bg  () { return "\033[0;40m"; }
    //! \brief Set the background color to red.
    inline std::string red_bg    () { return "\033[0;41m"; }
    //! \brief Set the background color to green.
    inline std::string green_bg  () { return "\033[0;42m"; }
    //! \brief Set the background color to yellow.
    inline std::string yellow_bg () { return "\033[0;43m"; }
    //! \brief Set the background color to blue.
    inline std::string blue_bg   () { return "\033[0;44m"; }
    //! \brief Set the background color to purple.
    inline std::string purple_bg () { return "\033[0;45m"   ; }
    //! \brief Set the background color to cyan.
    inline std::string cyan_bg   () { return "\033[0;46m"; }
    //! \brief Set the background color to white.
    inline std::string white_bg  () { return "\033[0;47m"; }

    //! \brief Resets the terminal, clearing any previous formatting command.
    inline std::string reset     () { return "\033[0m"; }
  }
}

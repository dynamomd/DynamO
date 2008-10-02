/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef IS_Stream_Op_H
#define IS_Stream_Op_H

#include <iostream>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>

#ifdef DYNAMO_Colour
# define IC_black "\033[22;30m" 
# define IC_red "\033[22;31m"
# define IC_green "\033[22;32m"
# define IC_blue "\033[22;34m"
# define IC_cyan "\033[22;36m"
# define IC_purple "\033[35m"
# define IC_white "\033[01;37m"
# define IC_white_brown "\033[43m\033[37m"
# define IC_blink "\033[5m"
# define IC_blink_off "\033[25m"
# define IC_reset "\033[0m"
#else
# define IC_black       ""
# define IC_red         ""
# define IC_green       ""
# define IC_blue        ""
# define IC_cyan        ""
# define IC_purple      ""
# define IC_white       ""
# define IC_white_brown ""
# define IC_blink       ""
# define IC_blink_off   ""
# define IC_reset       ""
#endif

#define IC_exception "\033[5m\033[41m\033[01;37m"

namespace DYNAMO
{
  inline std::string searchReplace(std::string in, const std::string& from, const std::string& to)
  {
    if (!in.empty())
      {
	std::string::size_type toLen = to.length();
	std::string::size_type frLen = from.length();
	std::string::size_type loc = 0;
	
	while (std::string::npos != (loc = in.find(from, loc)))
	  {
	    in.replace(loc, frLen, to);
	    loc += toLen;
	    
	    if (loc >= in.length())
	      break;
	  }
      }
    return in;
  }

  class Stream_Operator
  {
  public:
    Stream_Operator(const char* const & aName, const char* const & aColor):
      name(aName),color(aColor) 
    {
      OutputStream = &std::cout;
    };
    
    template<class T>
    const Stream_Operator& operator<<(T m) const
    { 
      std::ostringstream tmp("");
      tmp << m;
      *OutputStream << nReplace(tmp.str());
      return *this;
    }

    friend Stream_Operator operator <<(std::ostream &os, Stream_Operator SO)
    {
      SO.OutputStream = &os;
      return SO;
    }

    std::ostream& getStream() const
    { return *OutputStream; }
    
  protected:
    const char* const & name;
    const char* const & color;
    std::ostream *OutputStream;

    std::string nReplace(const std::string &message) const
    {
      return searchReplace(message, "\n", boost::lexical_cast<std::string>("\n") 
			   + boost::lexical_cast<std::string>(color) 
			   + boost::lexical_cast<std::string>(name) 
			   + boost::lexical_cast<std::string>(" :") 
			   + boost::lexical_cast<std::string>(IC_reset));
    }
    
  };

  class Colorise_Text_Stream_Operator
  {
  public:
    Colorise_Text_Stream_Operator(const char* const aColor):
      color(aColor)
    {
      OutputStream = &std::cout;
    };
    
    template<class T>
    const Colorise_Text_Stream_Operator& operator<<(T m) const
    { 
      *OutputStream << m;
      return *this;
    }
    
    const Colorise_Text_Stream_Operator& operator<<(const char *m) const
    { 
      *OutputStream << color << m << IC_reset;
      return *this;
    }


    friend Colorise_Text_Stream_Operator operator <<(std::ostream &os, Colorise_Text_Stream_Operator SO)
    {
      SO.OutputStream = &os;
      return SO;
    }

  protected:
    const char* const color;
    std::ostream *OutputStream;
  };
}

#endif

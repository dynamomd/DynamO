// (C) Copyright 2008 CodeRage, LLC (turkanis at coderage dot com)
// (C) Copyright 2005-2007 Jonathan Turkanis
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

// Adapted from an example of James Kanze. See
// http://www.gabi-soft.fr/codebase-en.html.

// See http://www.boost.org/libs/iostreams for documentation.

#ifndef BOOST_IOSTREAMS_LINE_WRAPPING_FILTER_HPP_INCLUDED
#define BOOST_IOSTREAMS_LINE_WRAPPING_FILTER_HPP_INCLUDED

#include <cstdio>                            // EOF.
#include <boost/iostreams/concepts.hpp>      // output_filter.
#include <boost/iostreams/operations.hpp>    // boost::iostreams::put.

namespace boost { namespace iostreams {
    class line_wrapping_output_filter : public output_filter {
    public:
      explicit line_wrapping_output_filter(int line_length = 80)
	: line_length_(line_length), col_no_(0)
      { }

      template<typename Sink>
      bool put(Sink& dest, int c)
      {
	if (c != '\n' && col_no_ >= line_length_ && !put_char(dest, '\n'))
	  return false;
	return put_char(dest, c);
      }

      template<typename Sink>
      void close(Sink&) { col_no_ = 0; }
    private:
      template<typename Sink>
      bool put_char(Sink& dest, int c)
      {
	if (!iostreams::put(dest, c))
	  return false;
	if (c != '\n')
	  ++col_no_;
	else
	  col_no_ = 0;
	return true;
      }
      int  line_length_;
      int  col_no_;
    };

  }
}

#endif

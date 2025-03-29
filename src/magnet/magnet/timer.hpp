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
#include <chrono>
#include <iostream>

namespace magnet {
/*! \brief A class to aid in timing the execution of commands.

  Example usages are

  \code
  {
    Timer timer("MyTestFunction()");
    for (size_t i(0); i < 1000; ++i, ++timer)
      MyTestFunction();
  }
  \endcode

  \code
  Timer timer();
  MyTestFunction();
  std::cout << "MyTestFunction() took " << timer.duration<std::micro>() << "
  microseconds";
  \endcode
*/
class Timer {
public:
  Timer(std::string text = "") : _text(text) {}

  ~Timer() {
    if (!_text.empty())
      std::cerr << _text << " "
                << duration<std::micro>() / (_count + (_count == 0))
                << " micro-s" << (_count ? " / call\n" : "\n");
  }

  void operator++() { ++_count; }

  template <class Period> double duration() const {
    return std::chrono::duration<double, Period>(
               std::chrono::steady_clock::now() - _start)
        .count();
  }

  std::chrono::time_point<std::chrono::steady_clock> _start =
      std::chrono::steady_clock::now();

private:
  size_t _count = 0;
  std::string _text;
};
} // namespace magnet

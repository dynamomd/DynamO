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
#include <magnet/exception.hpp>
#include <thread>
#include <vector>

namespace magnet {
namespace thread {
class ThreadGroup {
public:
  inline ThreadGroup() {}

  inline ~ThreadGroup() { join_all(); }

  template <class Function, class... Args>
  inline void create_thread(Function &&f, Args &&...args) {
    _threads.push_back(std::thread(f, args...));
  }

  inline size_t size() const { return _threads.size(); }

  inline void join_all() {
    for (auto &thread : _threads)
      if (thread.joinable())
        thread.join();
    _threads.clear();
  }

protected:
  std::vector<std::thread> _threads;

  ThreadGroup(const ThreadGroup &);
  ThreadGroup &operator=(const ThreadGroup &);
};
} // namespace thread
} // namespace magnet

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

#include <functional>
#include <mutex>
#include <queue>

namespace magnet {
namespace thread {
class TaskQueue {
public:
  // Actual queuer
  inline virtual void queueTask(std::function<void()> threadfunc) {
    std::lock_guard<std::mutex> lock(_queue_mutex);
    _waitingFunctors.push(threadfunc);
  }

  inline virtual void
  queueTasks(std::vector<std::function<void()>> &threadfuncs) {
    std::lock_guard<std::mutex> lock(_queue_mutex);

    for (const auto &func : threadfuncs)
      _waitingFunctors.push(func);

    threadfuncs.clear();
  }

  void drainQueue() {
    _queue_mutex.lock();

    while (!_waitingFunctors.empty()) {
      std::function<void()> task = _waitingFunctors.front();
      _waitingFunctors.pop();
      _queue_mutex.unlock();
      task();
      _queue_mutex.lock();
    }
    _queue_mutex.unlock();
  }

  virtual ~TaskQueue() {
    std::lock_guard<std::mutex> lock(_queue_mutex);
    _waitingFunctors = std::queue<std::function<void()>>();
  }

protected:
  std::queue<std::function<void()>> _waitingFunctors;
  std::mutex _queue_mutex;
};
} // namespace thread
} // namespace magnet

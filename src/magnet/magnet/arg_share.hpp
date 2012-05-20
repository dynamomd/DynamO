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
#include <cstring>
#include <magnet/exception.hpp>

namespace magnet {
  struct ArgShare {
    ArgShare():_argc(NULL), _argv(NULL) {}

    inline static ArgShare& getInstance() 
    {
      static ArgShare instance;
      return instance;
    }
    
    void setArgs(int& argc, char**& argv)
    {
      _argc = &argc;
      _argv = &argv;
    }

    int& getArgc()
    {
      if (_argc == NULL) M_throw() << "Command line args not passed to ArgShare";
      return *_argc;
    }

    char**& getArgv()
    {
      if (_argv == NULL) M_throw() << "Command line args not passed to ArgShare";

      return *_argv;
    }

  private:
    int* _argc;
    char*** _argv;
  };
}

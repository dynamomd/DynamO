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

#ifndef EventTypeTracking_H
#define EventTypeTracking_H

#include <utility>
#include <string>
#include "../dynamics/eventtypes.hpp"

namespace DYNAMO
{
  struct SimData;
}

class Interaction;
class GlobalEvent;
class System;
class LocalEvent;
class IntEvent;

namespace EventTypeTracking {
  
  //! Keeps the type of event (interaction etc) and the ID num
  typedef std::pair<size_t, EEventType> classKey;

  std::string getName(const classKey&, const DYNAMO::SimData*);

  classKey getClassKey(const IntEvent&);

  classKey getClassKey(const System&);

  classKey getClassKey(const GlobalEvent&);

  classKey getClassKey(const LocalEvent&);
}

#endif

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

#pragma once

#ifndef EventTypeTracking_H
#define EventTypeTracking_H

#include <utility>
#include <string>

namespace DYNAMO
{
  struct SimData;
}

class CInteraction;
class CGlobal;
class CSystem;

namespace EventTypeTracking {
  //! \brief This is to stop the use of maps
  enum eventClass
    {
      NOEventClass, //!< This is the initial last event type for particles
      InteractionClass,
      GlobalClass,
      SystemClass
    };
  
  typedef std::pair<size_t, eventClass> classKey;

  std::string getName(const classKey&, const DYNAMO::SimData*);

  classKey getClassKey(const CInteraction&);

  classKey getClassKey(const CSystem&);

  classKey getClassKey(const CGlobal&);
}

#endif

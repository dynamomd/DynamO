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
#include "eventtypetracking.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/include.hpp"

namespace EventTypeTracking {

  std::string getName(const classKey& key, const DYNAMO::SimData* Sim)
  {
  switch (key.second)
    {
    case INTERACTION:
      return Sim->dynamics.getInteractions()[key.first]->getName();
      break;
    case GLOBAL:
      return Sim->dynamics.getGlobals()[key.first]->getName();
      break;
    case SYSTEM:
      return Sim->dynamics.getSystemEvents()[key.first]->getName();
      break;
    case LOCAL:
      return Sim->dynamics.getLocals()[key.first]->getName();
      break;
    default:
      M_throw() << "Collision matrix found an unknown event class";
    }
  }

  classKey getClassKey(const IntEvent& i)
  {
    return classKey(i.getInteractionID(), INTERACTION);
  }
  classKey getClassKey(const System& s)
  {
    return classKey(s.getID(), SYSTEM);
  }

  classKey getClassKey(const GlobalEvent& g)
  {
    return classKey(g.getGlobalID(), GLOBAL);
  }

  classKey getClassKey(const LocalEvent& g)
  {
    return classKey(g.getLocalID(), LOCAL);
  }
}

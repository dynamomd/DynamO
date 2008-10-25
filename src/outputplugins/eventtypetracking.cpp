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
#include "eventtypetracking.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/include.hpp"

namespace EventTypeTracking {

  std::string getName(const classKey& key, const DYNAMO::SimData* Sim)
  {
  switch (key.second)
    {
    case InteractionClass:
      return Sim->Dynamics.getInteractions()[key.first]->getName();
      break;
    case GlobalClass:
      return Sim->Dynamics.getGlobals()[key.first]->getName();
      break;
    case SystemClass:
      return Sim->Dynamics.getSystemEvents()[key.first]->getName();
      break;
    default:
      D_throw() << "Collision matrix found an unknown event class";
    }
  }

  classKey getClassKey(const CInteraction& i)
  {
    return classKey(i.getID(), InteractionClass);
  }
  classKey getClassKey(const CSystem& s)
  {
    return classKey(s.getID(), SystemClass);
  }

  classKey getClassKey(const CGlobEvent& g)
  {
    return classKey(g.getGlobalID(), GlobalClass);
  }

  classKey getClassKey(const CLocalEvent& g)
  {
    return classKey(g.getLocalID(), GlobalClass);
  }
}

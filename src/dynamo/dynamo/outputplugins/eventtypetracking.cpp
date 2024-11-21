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
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/include.hpp>

namespace dynamo {
  namespace EventTypeTracking {

    std::string getEventSourceName(const EventSourceKey& key, const dynamo::Simulation* Sim)
    {
      switch (key.second)
	{
	case INTERACTION:
	  return Sim->interactions[key.first]->getName();
	  break;
	case GLOBAL:
	  return Sim->globals[key.first]->getName();
	  break;
	case SYSTEM:
	  return Sim->systems[key.first]->getName();
	  break;
	case LOCAL:
	  return Sim->locals[key.first]->getName();
	  break;
	default:
	  M_throw() << "Collision matrix found an unknown event class";
	}
    }

    std::string getEventSourceTypeName(const EventSourceKey& key)
    {
      switch (key.second)
	{
	case INTERACTION: return "Interaction";
	case GLOBAL: return "Global";
	case SYSTEM: return "System";
	case LOCAL: return "Local";
	default:
	  M_throw() << "Collision matrix found an unknown event class";
	}
    }

    EventSourceKey getEventSourceKey(const Event& i)
    {
      return EventSourceKey(i._sourceID, i._source);
    }
  }
}

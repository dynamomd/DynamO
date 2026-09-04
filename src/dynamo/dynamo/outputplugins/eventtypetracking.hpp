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
#include <dynamo/eventtypes.hpp>
#include <string>
#include <utility>

namespace dynamo {
class Simulation;
class Interaction;
class System;

namespace EventTypeTracking {

//! Keeps the type of the event source (GLOBAL, LOCAL, INTERACTION, SYSTEM) and the ID
typedef std::pair<size_t, EventSource> EventSourceKey;

//! Combines the EventSourceKey (type i.e. INTERACTION, and ID) with the EEventType (CORE, WELL, WALL, VIRTUAL, etc)
typedef std::pair<EventSourceKey, EEventType> EventKey;

std::string getEventSourceName(const EventSourceKey &,
                               const dynamo::Simulation *);
std::string getEventSourceTypeName(const EventSourceKey &);

EventSourceKey getEventSourceKey(const Event &);
} // namespace EventTypeTracking
} // namespace dynamo

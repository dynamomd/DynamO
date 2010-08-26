/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#define ETYPE_ENUM_FACTORY(F)\
  F(NONE, /*!< No collision occurs*/ ) \
  F(CELL, /*!< Marks cell transitions*/ ) \
  F(GLOBAL, /*!< Marks all global events*/ ) \
  F(INTERACTION, /*!< Marks Interaction events*/ ) \
  F(SYSTEM, /*!< Marks System events*/ ) \
  F(LOCAL, /*!< Marks Local events*/ ) \
  F(CORE, /*!< Hard core collision*/ ) \
  F(WELL_IN, /*!< Well Event where spheres are heading toward each other*/ ) \
  F(WELL_OUT, /*!< Well Event where spheres are heading away from each other*/ ) \
  F(WELL_KEUP, /*!< Well event where Kinetic Energy increases*/ ) \
  F(WELL_KEDOWN, /*!< Well event where Kinetic Energy decreases*/ ) \
  F(BOUNCE, /*!< CORE event due to energetic constraints*/ ) \
  F(WALL, /*!< Wall or other obstacle event*/ ) \
  F(GAUSSIAN, /*!< Reassignment from a gaussian Andersen thermostat*/ ) \
  F(DSMC, /*!< DSMC event*/ ) \
  F(UMBRELLA, /*!< Umbrella potential event*/ ) \
  F(HALT, /*!< Call to halt the system*/ ) \
  F(STREAM, /*!< Call to free stream the system an amount*/ ) \
  F(NON_EVENT, /*!< Anything that is not part of the system dynamics*/ ) \
  F(RESCALE, /*!< A rescaling of the system energy*/ ) \
  F(VIRTUAL, /*!< Fake events that cause a particle to free stream*/ )

#define buildEnum(VAL,COMMENT) \
  VAL, COMMENT

typedef enum {
ETYPE_ENUM_FACTORY(buildEnum)
} EEventType; 

#include <ostream>

#define printEnum(VAL,COMMENT) \
  case VAL: return os << #VAL;

inline std::ostream& operator<<(std::ostream& os, EEventType etype)
{
  switch (etype)
    {
ETYPE_ENUM_FACTORY(printEnum)
    }

  return os << "Failed to find a name for the Event Type! Memory corruption?";
}

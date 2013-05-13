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
#include <ostream>
#include <magnet/exception.hpp>

namespace dynamo {
#define ETYPE_ENUM_FACTORY(F)						\
  F(NONE) /*!< No collision occurs*/ 					\
  F(CELL) /*!< Marks cell transitions */ 				\
  F(GLOBAL)  /*!< Marks all global events*/ 				\
  F(INTERACTION) /*!< Marks Interaction events*/ 			\
  F(SYSTEM) /*!< Marks System events*/ 				\
  F(LOCAL)  /*!< Marks Local events*/ 				\
  F(CORE) /*!< Hard core collision*/					\
  F(STEP_IN) /*!< Event where particles are heading toward each other and moved in a step of the potential*/ \
  F(STEP_OUT) /*!< Event where particles are heading away from each other and moved out a step of the potential*/ \
  F(NBHOOD_IN)  /*!< Event where particles enter a neighbourhood*/ 	\
  F(NBHOOD_OUT) /*!< Event where particles leave a neighbourhood*/	\
  F(BOUNCE) /*!< CORE event due to energetic constraints*/		\
  F(WALL) /*!< Wall or other obstacle event*/				\
  F(GAUSSIAN) /*!< Reassignment from a gaussian Andersen thermostat*/	\
  F(DSMC) /*!< DSMC event*/						\
  F(UMBRELLA) /*!< Umbrella potential event*/				\
  F(NON_EVENT) /*!< Anything that is not part of the system dynamics*/ \
  F(RESCALE) /*!< A rescaling of the system energy*/  		\
  F(RECALCULATE) /*!< Fake events that cause a particle to free stream*/ \
  F(RECALCULATE_PARABOLA) /*!< Fake event used to track when a particle goes through its parabola. Needed to keep the dynamics deterministic.*/ \
  F(VIRTUAL) /*< Passed to output plugins to let them know that this event is not a true event.*/ \
  F(ROTATEGRAVITY) /*!< Event to rotate the gravity vector*/ \
  F(SLEEP) /*!< Event to transition a particle from dynamic to static*/ \
  F(RESLEEP) /*!< Event to zero a sleeping particles velocity after being hit*/ \
  F(WAKEUP) /*!< Event to transition a particle from static to dynamic*/ \
  F(CORRECT) /*!< An event used to correct a previous event*/
  
#define buildEnum(VAL) VAL,
#define printEnum(VAL) case VAL: return os << #VAL;

  typedef enum {
    ETYPE_ENUM_FACTORY(buildEnum)
    FINAL_ENUM_TO_CATCH_THE_COMMA
  } EEventType; 

  inline std::ostream& operator<<(std::ostream& os, EEventType etype)
  {
    switch (etype)
      {
	ETYPE_ENUM_FACTORY(printEnum)
      case FINAL_ENUM_TO_CATCH_THE_COMMA:
      default:
	break;
      }

    M_throw() << "Failed to find a name for the Event Type! The value must be uninitialised or memory is being corrupted somehow.";
  }

#undef ETYPE_ENUM_FACTORY
#undef buildEnum
#undef printEnum
}


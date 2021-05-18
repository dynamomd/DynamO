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
#include <magnet/exception.hpp>
#include <algorithm>
#include <ostream>

namespace dynamo {
#define ETYPE_ENUM_FACTORY(F)						\
  F(NONE) /*!< No collision occurs*/ 					\
  F(CELL) /*!< Marks cell transitions */ 				\
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
  F(NON_EVENT) /*!< Anything that is not part of the system dynamics*/	\
  F(RESCALE) /*!< A rescaling of the system energy*/			\
  F(RECALCULATE) /*!< Fake events that cause a particle to free stream*/ \
  F(RECALCULATE_PARABOLA) /*!< Fake event used to track when a particle goes through its parabola. Needed to keep the dynamics deterministic.*/ \
  F(VIRTUAL) /*< Passed to output plugins to let them know that this event is not a true event.*/ \
  F(ROTATEGRAVITY) /*!< Event to rotate the gravity vector*/		\
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
    
  typedef enum {
    INTERACTION,
    LOCAL,
    GLOBAL,
    SYSTEM,
    SCHEDULER,
    NOSOURCE
  } EventSource; 

  inline std::ostream& operator<<(std::ostream& os, EventSource etype)
  {
    switch (etype)
      {
      case INTERACTION: os << "Interaction"; break;
      case LOCAL: os << "Local"; break;
      case GLOBAL: os << "Global"; break;
      case SYSTEM: os << "System"; break;
      case SCHEDULER: os << "Scheduler"; break;
      case NOSOURCE: os << "No-Source"; break;
      default:
	M_throw() << "Failed to find a name for the Event source type! The value must be uninitialised or memory is being corrupted somehow.";
	break;
      };
    return os;
  }

  /*! \brief A generic event type, which the more specialised events
    are converted to before they are sorted.

    This conversion is lossy, so events need to be recalculated if
    they are to be exectuted. 

    The RECALCULATE event type is special. If any IntEvent, GlobalEvent
    or LocalEvent has a type RECALCULATE, it is carried through. RECALCULATE
    events cause the system to be moved forward in time and the
    events for the particle are recalculated. This can all be
    handled by the scheduler.
  */
  class Event
  {
  public:   
    double _dt;
    size_t _particle1ID;
    size_t _sourceID;

    union {
      size_t _particle2ID;
      size_t _additionalData1;
    };

     union {
      size_t _particle2eventcounter;
      size_t _additionalData2;
    };
      
    EventSource _source;
    EEventType _type;
      
    inline Event():
      _dt(std::numeric_limits<float>::infinity()),
      _particle1ID(std::numeric_limits<size_t>::max()),
      _sourceID(std::numeric_limits<size_t>::max()),
      _additionalData1(std::numeric_limits<size_t>::max()),
      _additionalData2(std::numeric_limits<size_t>::max()),
      _source(NOSOURCE),
      _type(NONE)
    {}
      
      inline Event(size_t particle1ID, double dt, EventSource source, EEventType type, size_t sourceID, size_t additionalData1 = std::numeric_limits<size_t>::max(), size_t additionalData2 = std::numeric_limits<size_t>::max()) throw():
      _dt(dt),
      _particle1ID(particle1ID),
      _sourceID(sourceID),
      _additionalData1(additionalData1),
      _additionalData2(additionalData2),
      _source(source),
      _type(type)
    {}
    
    inline bool operator< (const Event& o) const throw()
    { return _dt < o._dt; }

    inline bool operator> (const Event& o) const throw()
    { return _dt > o._dt; }

    inline bool operator==(const Event& o) const throw()
    { 
      return (_dt == o._dt) && (_particle1ID == o._particle1ID) && 
	(_sourceID == o._sourceID) && (_type == o._type) && 
	(_additionalData1 == o._additionalData1) && (_source == o._source) && 
	((_source == INTERACTION) || (_additionalData2 == o._additionalData2));
    }
  };


  inline std::ostream& operator<<(std::ostream& os, Event event)
  {
    os << "Event{dt = " << event._dt << ", p1ID = " << event._particle1ID
       << ", sourceID = " << event._sourceID 
       << ", data1 = " << event._additionalData1
       << ", data2 = " << event._additionalData2
       << ", source = " << event._source
       << ", type = " << event._type
       << "}";
    return os;
  }
}


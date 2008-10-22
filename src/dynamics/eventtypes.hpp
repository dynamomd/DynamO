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

#ifndef CollTypes_H
#define CollTypes_H

//! \brief Event that occured
typedef enum {
  /* These should not be seen outside the scheduler*/
  CELL        , /*!< Marks cell transitions in the scheduler*/
  GLOBAL      , /*!< Marks all global events in the scheduler*/
  INTERACTION , /*!< Marks Interaction events in the scheduler*/
  /* These are real types of events outside the scheduler */
  NONE        , /*!< No collision occurs*/
  CORE        , /*!< Hard core collision*/
  WELL_IN     , /*!< Well Event, could be BOUNCE WELL_KEUP or WELL_KEDOWN*/
  WELL_OUT    , /*!< Well Event, could be BOUNCE WELL_KEUP or WELL_KEDOWN*/
  WELL_KEUP   , /*!< Well energy change where KE increases*/
  WELL_KEDOWN , /*!< Well energy change where KE decreases*/
  BOUNCE      , /*!< failed to exit/enter a well*/
  WALL        , /*!< Wall event*/
  GAUSSIAN    , /*!< Reassignment from a gaussian*/
  HALT        , /*!< Call to halt the system*/
  STREAM      , /*!< Call to free stream the system an amount*/
  NON_EVENT,    /*!< Anything like a ticker, that is not part of the system dynamics*/
  VIRTUAL       /*!< This is not an event yet it requires a
                   recalculation of the particles collision list*/
} EEventType;

#endif

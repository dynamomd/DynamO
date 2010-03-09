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

#ifndef CollTypes_H
#define CollTypes_H

//! \brief Event that occured
typedef enum {
  /* Generic names for events occuring*/
  NONE        , /*!< No collision occurs*/
  CELL        , /*!< Marks cell transitions.*/
  GLOBAL      , /*!< Marks all global events.*/
  INTERACTION , /*!< Marks Interaction events.*/
  SYSTEM      , /*!< Marks System events.*/
  LOCAL       , /*!< Marks Local events.*/
  /* More specific flags for event types. */
  CORE        , /*!< Hard core collision. */
  WELL_IN     , /*!< Well Event, while spheres are heading toward each
		  other. Could be BOUNCE WELL_KEUP or WELL_KEDOWN.*/
  WELL_OUT    , /*!< Well Event, while spheres are heading away from
		  each other could be BOUNCE WELL_KEUP or
		  WELL_KEDOWN. */
  WELL_KEUP   , /*!< Well event where Kinetic Energy increases. */
  WELL_KEDOWN , /*!< Well event where Kinetic Energy decreases. */
  BOUNCE      , /*!< Well event where the particles failed to change
                   their kinetic energy due to some constraint. */
  WALL        , /*!< Wall or other obstacle event. */
  GAUSSIAN    , /*!< Reassignment from a gaussian, Andersen thermostat. */
  DSMC        , /*!< DSMC event*/
  UMBRELLA    , /*!< Umbrella potential event*/
  HALT        , /*!< Call to halt the system. */
  STREAM      , /*!< Call to free stream the system an amount. */
  NON_EVENT,    /*!< Anything like a ticker, that is not part of the
		  system dynamics. Does not require an update of the
		  system in any way. */
  VIRTUAL       /*!< This is not an event yet it requires a
                   recalculation of the particles collision
                   list. Possibly used in a sentinal */ 
} EEventType;

#endif

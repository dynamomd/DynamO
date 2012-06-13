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
/*! \file simulation.hpp
 * \brief Contains the definition of the Simulation class.
 */
#pragma once

#include <dynamo/base.hpp>
#include <dynamo/simdata.hpp>

namespace dynamo {
  class Dynamics;
  class OutputPlugin;

  /*! \brief A single simulation holding particles, dynamics and output
   * plugins.
   *
   * This class is the typical realisation of a simulation program. It
   * can pretty much perform a standard simulation without any other
   * supporting class structure like the Engine and Coordinator. This
   * class handles the interface to the simulation and also stores the
   * Simulation data by deriving from the dynamo::SimData class.
   *
   *
   */
  class Simulation: public dynamo::SimData
  {
  public:
  };
}

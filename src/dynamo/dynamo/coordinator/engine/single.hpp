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
/*! \file single.hpp
 * \brief Contains the simulation Engine ESingleSimulation
 */

#pragma once
#include <dynamo/coordinator/engine/engine.hpp>

namespace dynamo {
  /*! \brief An Engine for simulating a single system.
   *
   * This merely sets up and exectues a single Simulation instance.
   */
  class ESingleSimulation: public Engine
  {
  public:
    /*! \brief Only constructor.
     *
     * \param vm A reference to the Coordinator's parsed command line variables.
     * \param tp A reference to the thread pool of the dynarun instance.
     */ 
    ESingleSimulation(const boost::program_options::variables_map& vm, 
		      magnet::thread::ThreadPool& tp);

    /*! \brief Trivial virtual destructor */
    virtual ~ESingleSimulation() {}
  
    /*! \brief There is no status to be printed other than what the
     * Simulation outputs.
     */
    virtual void printStatus() {}
  
    /*! \brief Just runs the Simulation::runSimulation() loop and
     * provides a peek functionality.
     */
    virtual void runSimulation();

    /*! \brief Just wraps the Simulation::outputData(const char*) function.
     */
    virtual void outputData();
  
    /*! \brief Just wraps the Simulation::writeXMLfile(const char*) function.
     */
    virtual void outputConfigs();

    /*! \brief No Engine finalisation required.
     */
    virtual void finaliseRun() {}

    /*! \brief Just wraps the Simulation::simShutdown() function.
     */
    virtual void forceShutdown();

    /*! \brief Performs the minimum steps to initialise a simulation.
     */
    virtual void initialisation();

    /*! \brief Triggers the peek mode in the loop.
     */
    virtual void peekData();

  protected:
    /*! \brief The single instance of a Simulation required.
     */
    Simulation simulation;
  
    /*! \brief If this is true, the Simulation's end time will be reset
     * and the ESingleSimulation run loop in runSimulation() will be repeated.
     */
    bool peekMode; 
  };
}

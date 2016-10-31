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

/*! \file coordinator.hpp
 *
 * \brief Contains the header code for the Coordinator class.
 */

#pragma once

#include <dynamo/coordinator/engine/engine.hpp>
#include <magnet/thread/threadpool.hpp>
#include <boost/program_options.hpp>
#include <vector>

namespace dynamo {

  /*! \brief The main class for the dynarun program.
   
    This class is responsible for sorting out the correct simulation Engine to 
    run and initialising computational node specific objects like the 
    ThreadPool.
   */
  class Coordinator
  {
  private:

    /*! \brief The default constructor.

      This constructor is hidden as part of the singleton nature of
      this coordinator (i.e., there can only ever be one coordinator
      in a single program).
     */
    Coordinator() {}

  public:
    void enableVisualisation() { _enableVisualisation = true; }


    /*! \brief This is how the singleton Coordinator class is
        accessed.
     */
    static Coordinator& get() 
    { 
      static Coordinator _coordinator;
      return _coordinator; 
    }

    /*! \brief Parses the command line options, including any engine specific options.
     
      This function must know how to get the command line options for
      all available engines. 
     
      \return The parsed options vm is returned in case the owning
      class/function needs to inspect them.
     
      \param argc The number of command line arguments.
      \param argv A pointer to the array of command line arguments.
     */
    boost::program_options::variables_map& parseOptions(int argc, char* argv[]);

    /*! \brief Creates the specified Engine according to the command
      line options and initialises it.
     */
    void initialise();

    /*! \brief Calls Engine::runSimulation() if there are collisions to execute.
     */
    void runSimulation();

    /*! \brief Outputs any simulation data collected using Engine::outputData().
     
      In the future this will also output any data collected on the engine or 
      system state, i.e. the mpi subsystem.
     */ 
    void outputData();
  
    /*! \brief Calls Engine::outputConfigs() to print the final configurations
     * if any dynamics was actually run.
     */
    void outputConfigs();

    /*! \brief The signal handler for the dynarun.
     
      The purpose of this function is to respond to singals
      gracefully. If the user presses Ctrl-c they will be presented with a menu
      describing the options available.
     
      When the program is run in a batch control system like PBS or SGE
      and the job approaches its time limits the queuing system sends
      SIGUSR1 and SIGUSR2 shortly before to allow the program to
      gracefully exit. We catch these signals and shutdown as quickly
      as possible.
     */
    static void setup_signal_handler();

    /*! \brief The signal handler for the dynarun.
     
      The purpose of this function is to respond to singals
      gracefully. If the user presses Ctrl-c they will be presented with a menu
      describing the options available.
     
      When the program is run in a batch control system like PBS or SGE
      and the job approaches its time limits the queuing system sends
      SIGUSR1 and SIGUSR2 shortly before to allow the program to
      gracefully exit. We catch these signals and shutdown as quickly
      as possible.
     */
#ifdef _WIN32
    static BOOL signal_handler(DWORD);
#else
    static void signal_handler(int);
#endif
    
  private:
    /*! \brief Contains the parsed command line options, engines carry references
      to these values.
     */
    boost::program_options::variables_map vm;

    /*! \brief A smart pointer to the Engine being run.
     */
    shared_ptr<Engine> _engine;

    /*! \brief A thread pool to utilise multiple cores on the
      computational node.
      
      This ThreadPool is used/referenced by all code in a single
      dynarun process.
    */
    magnet::thread::ThreadPool _threads;
    
    bool _enableVisualisation;
  };
}

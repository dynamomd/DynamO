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
/*! \file engine.hpp
 * \brief Contains the definition of Engine.
 */

#pragma once

#include <boost/program_options.hpp>
#include <dynamo/simulation.hpp>

namespace magnet {
namespace thread {
class ThreadPool;
}
} // namespace magnet

namespace dynamo {
/*! \brief An engine to control/manipulate one or more Simulation's.
 *
 * Engine is a virtual base class interface for many different
 * engines. These engines manipulate Simulation(s) data by running
 * them and/or altering them for the purpose of a study.
 *
 * The simplest engine is ESingleSimulation and probably the best one to try
 * and understand at first.
 *
 * The initialisation() steps of an engine have been broken up into
 * three stages so that the derived engines can hook in where they
 * need to.
 * - preSimInit() - Before the Simulation(s) is/are initialised
 * - setupSim(Simulation &, const std::string) - The initialisation of the
 * Simulation(s)
 * - postSimInit() - After the Simulation(s) is/are initialised
 *
 */
class Engine {
public:
  /*! \brief The default constructor.
   *
   * \param vm Reference to the parsed command line variables.
   * \param configFile A format string on how config files should be written
   * out.
   * \param outputFile A format string on how output files should be written
   * out.
   * \param tp The processes ThreadPool for parallel processing.
   */
  Engine(const boost::program_options::variables_map &vm,
         std::string configFile, std::string outputFile,
         magnet::thread::ThreadPool &tp);

  /*! \brief The trivial virtual destructor. */
  virtual ~Engine() {}

  /*! \brief A hook for the initialisation stage of an engine
   *
   * This function should at the very least call in the following order
   * - preSimInit()
   * - setupSim(Simulation &, const std::string) for every Simulation class
   * - postSimInit()
   *
   */
  virtual void initialisation() = 0;

  /*! \brief This hook is run before the engine is destroyed.
   *
   * This is if the engine needs to change its state before shutting
   * down.
   * E.g. the ECompressingSimulation needs to change the Dynamics
   * back to the old one.
   */
  virtual void finaliseRun() = 0;

  /*! \brief The main simulation "loop"/call for the engine
   *
   * Some engines like the EReplicaExchangeSimulation require a loop and it will
   * be implemented here
   */
  virtual void runSimulation() = 0;

  /*! \brief Output any data collected during the run by the
   * Simulation's and the Engine.
   */
  virtual void outputData() = 0;

  /*! \brief Output the configurations of the Simulation's and Engine
   * so the run can be continued.
   *
   * This function must be safe to call during an interrupt.
   */
  virtual void outputConfigs() = 0;

  /*! \brief Add common options for all the engines to the options_description
   *
   * Each engine will define a similar static function to add their
   * options.
   *
   * \param od The boost options description that the common engine
   * options are to be added to.
   */
  static void getCommonOptions(boost::program_options::options_description &od);

  void sigint() { _SIGINT = true; }
  void sigterm() { _SIGTERM = true; }

  void enableVisualisation() { _loadVisualiser = true; }

protected:
  /*! \brief Code common to most engines pre simulation initialisation.
   */
  virtual void preSimInit();

  /*! \brief Code common to loading a Simulation from a config file.
   *
   * \param Sim Simulation to set up.
   * \param inFile Name of configuration file to load.
   */
  virtual void setupSim(Simulation &Sim, const std::string inFile);

  /*! \brief Once the Simulation is loaded and initialised you may
   * need to alter it/load plugins/initialise some Engine datastruct.
   */
  virtual void postSimInit(Simulation &) {}

  /*! \brief A reference to the Coordinators parsed command line variables.
   */
  const boost::program_options::variables_map &vm;

  std::string configFormat;
  std::string outputFormat;
  bool _SIGINT;
  bool _SIGTERM;
  bool _loadVisualiser;
  magnet::thread::ThreadPool &threads;
};
} // namespace dynamo

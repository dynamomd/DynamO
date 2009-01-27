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
/*! \file engine.hpp
 * \brief Contains the definition of CEngine.
 */

#ifndef CEngine_H
#define CEngine_H

#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>
#include "../../simulation/simulation.hpp"

class CThreadPool;

/*! \brief An engine to control/manipulate one or more CSimulation's.
 *
 * CEngine is a virtual base class interface for many different
 * engines. These engines manipulate CSimulation(s) data by running
 * them and/or altering them for the purpose of a study.
 * 
 * The simplest engine is CESingle and probably the best one to try
 * and understand at first.
 *
 * The initialisation() steps of an engine have been broken up into
 * three stages so that the derived engines can hook in where they
 * need to.
 * - preSimInit() - Before the CSimulation(s) is/are initialised
 * - setupSim(CSimulation &, const std::string) - The initialisation of the CSimulation(s)
 * - postSimInit() - After the CSimulation(s) is/are initialised
 *
 */
class CEngine
{
public:
  /*! \brief The default constructor.
   *
   * \param vm Reference to the parsed command line variables.
   * \param configFile A format string on how config files should be written out.
   * \param outputFile A format string on how output files should be written out.
   * \param tp The processes CThreadPool for parallel processing.
   */
  CEngine(const boost::program_options::variables_map& vm,
	  std::string configFile, std::string outputFile,
	  CThreadPool& tp);
  
  /*! \brief The trivial virtual destructor. */
  virtual ~CEngine() {}

  /*! \brief A hook for the initialisation stage of an engine
   * 
   * This function should at the very least call in the following order
   * - preSimInit()
   * - setupSim(CSimulation &, const std::string) for every CSimulation class
   * - postSimInit()
   *
   */
  virtual void initialisation() = 0;

  /*! \brief This hook is run before the engine is destroyed.
   *
   * This is if the engine needs to change its state before shutting
   * down.  
   * E.g. the CECompressor needs to change the CLiouvillean
   * back to the old one.
   */
  virtual void finaliseRun() = 0;

  /*! \brief Try to shut the engine down prematurely due to an
   * interrupt being called.
   *
   * This function must be safe to call during an interrupt.
   */
  virtual void forceShutdown() = 0;

  /*! \brief Called when the user requests the status of the currently
   * running engine.
   *
   * Output minimal data that you would like to track here.
   * This function must be safe to call during an interrupt.
   */
  virtual void printStatus() = 0;

  /*! \brief The main simulation "loop"/call for the engine
   *
   * Some engines like the CEReplexer require a loop and it will be
   * implemented here
   */
  virtual void runSimulation() = 0;

  /*! \brief Output any data collected during the run by the
   * CSimulation's and the CEngine.
   */
  virtual void outputData() = 0;

  /*! \brief Instruct the system to output its data using outputData()
   * at the next available point, for mid simulation previews.
   */
  virtual void peekData() = 0;

  /*! \brief Output the configurations of the CSimulation's and CEngine
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
  static void getCommonOptions(boost::program_options::options_description& od);
  
protected:
  /*! \brief Code common to most engines pre simulation initialisation.
   */
  virtual void preSimInit();

  /*! \brief Code common to loading a CSimulation from a config file.
   *
   * \param Sim CSimulation to set up.
   * \param inFile Name of configuration file to load.
   */
  virtual void setupSim(CSimulation & Sim, const std::string inFile);

  /*! \brief Once the CSimulation is loaded and initialised you may
   * need to alter it/load plugins/initialise some CEngine datastruct.
   */
  virtual void postSimInit(CSimulation&) {}

  /*! \brief A reference to the CCoordinators parsed command line variables.
   */
  const boost::program_options::variables_map& vm;
    
  std::string configFormat;
  std::string outputFormat;

  CThreadPool& threads;
};

#endif

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
#include "../../extcode/threadpool.hpp"

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
   */
  CEngine(const boost::program_options::variables_map& vm,
	  std::string configFile, std::string outputFile);
  
  /* \brief The trivial virtual destructor. */
  virtual ~CEngine() {}

  /* \brief A hook for the initialisation stage of an engine
   * 
   * This function should at the very least call in the following order
   * - preSimInit()
   * - setupSim(CSimulation &, const std::string) for every CSimulation class
   * - postSimInit()
   *
   */
  virtual void initialisation() = 0;

  /* \brief This hook is run before the engine is destroyed.
   *
   * If the engine has data to output it should do it in this function.
   */
  virtual void finaliseRun() = 0;

  /* \brief Try to shut the engine down prematurely.
   */
  virtual void forceShutdown() = 0;

  /* \brief Called when the user requests the status of the currently running engine.
   *
   * Output minimal data that you would like to track here
   */
  virtual void printStatus() = 0;

  virtual void runSimulation() = 0;

  virtual void outputData() = 0;

  virtual void outputConfigs() = 0;
  
  virtual void peekData() = 0;
  
  static void getCommonOptions(boost::program_options::options_description&);
  
protected:
  virtual void preSimInit();

  virtual void setupSim(CSimulation &, const std::string);

  virtual void postSimInit(CSimulation&);

  const boost::program_options::variables_map& vm;
  CThreadPool threads;
  std::string configFormat;
  std::string outputFormat;
};

#endif

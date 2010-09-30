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
/*! \file compressor.hpp
 * Contains the definition of ECompressingSimulation.
 */
#ifndef ECompressingSimulation_H
#define ECompressingSimulation_H

#include "single.hpp"
#include "../../inputplugins/compression.hpp"

/*! \brief This Engine compresses a configuration using the
 * LCompression Liouvillean.
 *
 * This is essentially a ESingleSimulation but with some extra steps to load
 * the compression Liouvillean at the start and then to restore the
 * old Liouvillean at the end.
 */
class ECompressingSimulation: public ESingleSimulation
{
public:
  /*!\brief The only constructor.
   *
   * \param vm The parsed command line options.
   * \param tp The shared thread pool.
   */
  ECompressingSimulation(const boost::program_options::variables_map& vm,
			 magnet::thread::ThreadPool& tp);

  /*! \brief A trivial virtual destructor
   */
  virtual ~ECompressingSimulation() {}

  /*! \brief Load the original Liouvillean before outputing the
   * configurations.
   * 
   * This is one of the few classes that does need to finalise before
   * output to restore the original system at a higher density.
   */
  virtual void finaliseRun();

  /*! \brief The options specific to the ECompressingSimulation class.
   *
   * This is used by the Coordinator::parseOptions function.
   *
   * \param od The options description to add the ECompressingSimulation options to.
   */
  static void getOptions(boost::program_options::options_description& od);

protected:
  /*! \brief Boot a CIPCompression plugin to handle the manipulation
   * of the single Simulation.
   *
   * This function also calls the Engine::preSimInit function
   */
  virtual void preSimInit();
  
  /*! \brief Use the CIPCompression plugins to switch to compression
   * dynamics.
   */
  virtual void setupSim(Simulation&, const std::string);

  /*! \brief A single CIPCompression plugin to manipulate the Simulation.
   */
  ClonePtr<CIPCompression> compressPlug;
};

#endif

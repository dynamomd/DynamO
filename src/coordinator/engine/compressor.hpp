/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
 * Contains the definition of CECompressor.
 */
#ifndef CECompressor_H
#define CECompressor_H

#include "single.hpp"
#include "../../inputplugins/compression.hpp"

/*! \brief This CEngine compresses a configuration using the
 * LCompression Liouvillean.
 *
 * This is essentially a CESingle but with some extra steps to load
 * the compression Liouvillean at the start and then to restore the
 * old Liouvillean at the end.
 */
class CECompressor: public CESingle
{
public:
  /*!\brief The only constructor.
   *
   * \param vm The parsed command line options.
   * \param tp The shared thread pool.
   */
  CECompressor(const boost::program_options::variables_map& vm,
	       CThreadPool& tp);

  /*! \brief A trivial virtual destructor
   */
  virtual ~CECompressor() {}

  /*! \brief Load the original Liouvillean before outputing the
   * configurations.
   * 
   * This is one of the few classes that does need to finalise before
   * output to restore the original system at a higher density.
   */
  virtual void finaliseRun();

  /*! \brief The options specific to the CECompressor class.
   *
   * This is used by the CCoordinator::parseOptions function.
   *
   * \param od The options description to add the CECompressor options to.
   */
  static void getOptions(boost::program_options::options_description& od);

protected:
  /*! \brief Boot a CIPCompression plugin to handle the manipulation
   * of the single CSimulation.
   *
   * This function also calls the CEngine::preSimInit function
   */
  virtual void preSimInit();
  
  /*! \brief Use the CIPCompression plugins to switch to compression
   * dynamics.
   */
  virtual void setupSim(CSimulation&, const std::string);

  /*! \brief A single CIPCompression plugin to manipulate the CSimulation.
   */
  smrtPlugPtr<CIPCompression> compressPlug;
};

#endif

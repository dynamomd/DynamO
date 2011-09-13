/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
/*! \file compression.hpp
 * Contains the definition of the class CIPCompression.
 */

#pragma once
#include <dynamo/inputplugins/inputplugin.hpp>

namespace dynamo {
  class Liouvillean;

  /*! \brief A plugin to change a simulation to compression dynamics and
   * back again.
   * 
   * This class came about as when a simulation is being compressed its
   * dynamics, or more specifically its Liouvillean, is replaced with
   * the LCompression liouvilean. This stores the old liouvillean and
   * also provides several helpful plugins to hack parts of the system
   * into co-operating with the compression like the cellular scheduler.
   */
  class CIPCompression: public CInputPlugin
  {
  public:
  
    /*! \brief The only constructor.
     *
     * \param sim The Simulation this plugin is in control of
     * \param cr The compression rate of the Simulation.
     */
    CIPCompression(dynamo::SimData* sim, double cr);

    /*! \brief Stores the old CLiovillean and installs the LCompression.
     */  
    void MakeGrowth();
  
    /*! \brief Restores the old Liouvillean stored in oldLio.
     */
    void RestoreSystem();  

    /*! \brief Installs the CSCellHack system event to make sure the
     *   cellular scheduler doesn't fail.
     */
    void CellSchedulerHack();

    /*! \brief Limits the maximum packing fraction by installing a
     *  CStHalt system event at the right time.
     *
     * \param mxpf The maximum packing fraction allowed.
     */
    void limitPackingFraction(double mxpf);

    /*! \brief Limits the maximum density by installing a
     *  CStHalt system event at the right time.
     *
     * \param mxrho The maximum number density allowed.
     */
    void limitDensity(double mxrho);

  private:
    /*! \brief The compression rate of the simulation.
     */
    double growthRate;
  
    /*! \brief The old Liouvillean of the simulation.
     */
    std::tr1::shared_ptr<Liouvillean> oldLio;
  
    /*! \brief Stores a cell overlap parameter of the cellular scheduler
     * to be restored later.
     */
  };
}

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
/*! \file compression.hpp
 * Contains the definition of the class CIPCompression.
 */

#ifndef CIPCompression_H
#define CIPCompression_H

#include "inputplugin.hpp"

class CLiouvillean;

/*! \brief A plugin to change a simulation to compression dynamics and
 * back again.
 * 
 * This class came about as when a simulation is being compressed its
 * dynamics, or more specifically its CLiouvillean, is replaced with
 * the CLCompression liouvilean. This stores the old liouvillean and
 * also provides several helpful plugins to hack parts of the system
 * into co-operating with the compression like the cellular scheduler.
 */
class CIPCompression: public CInputPlugin
{
 public:
  
  /*! \brief The only constructor.
   *
   * \param sim The CSimulation this plugin is in control of
   * \param cr The compression rate of the CSimulation.
   */
  CIPCompression(DYNAMO::SimData* sim, Iflt cr);

  /*! \brief Stores the old CLiovillean and installs the CLCompression.
   */  
  void MakeGrowth();
  
  /*! \brief Restores the old CLiouvillean stored in oldLio.
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
  void limitPackingFraction(Iflt mxpf);

  /*! \brief An expensive sanity check for the system.
   *
   * Ensures that the compression dynamics haven't corrupted the
   * system configuration by introducing invalid states.
   */
  void checkOverlaps();

private:
  /*! \brief The compression rate of the simulation.
   */
  Iflt growthRate;
  
  /*! \brief The old CLiouvillean of the simulation.
   */
  CLiouvillean* oldLio;
  
  /*! \brief Stores a cell overlap parameter of the cellular scheduler
   * to be restored later.
   */
  Iflt oldLambda;
};

#endif

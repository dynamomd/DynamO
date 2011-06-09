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
/*! \file simulation.hpp
 * \brief Contains the definition of the Simulation class.
 */
#pragma once

#include <dynamo/base.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <boost/scoped_array.hpp>

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
  /*! \brief Initialise the entire Simulation and the SimData struct.
   *
   * Most classes will have an initialisation function and its up to
   * this function to call them all and in the right order.
   */
  void initialise();

  //! Main loop for the Simulation
  //! \param silentMode If true, the periodic output of the simulation is supressed. 
  void runSimulation(bool silentMode = false);
  
  //! Writes the results of the Simulation to a file at the passed path.
  //! \param filename The path to the XML file to write (this file
  //! will either be created or overwritten). The filename must end in
  //! either ".xml" for uncompressed xml files or ".bz2" for bzip2
  //! compressed configuration files.
  void outputData(std::string filename = "output.xml.bz2");

  //! Used by dynamod to inform the Simulation that the configuration
  //! of the simulation is complete. This updates the Simulation
  //! status.
  //! \sa getStatus()
  void configLoaded();

  //! This function makes the Simulation exit the runSimulation loop
  //! at the next opportunity. Used when forcing the simulation to
  //! stop.
  void simShutdown();

  //! Sets how many events the Simulation loop should run for.
  //! \sa runSimulation
  void setTrajectoryLength(unsigned long long);

  //! Sets how many events should occur between outputting the
  //! Simulation state on the screen.
  //! \sa runSimulation
  void setnPrint(unsigned long long);

  //! Sets the random seed used by the Simulation random number
  //! generator.
  //! \sa dynamo::SimData::ranGenerator
  void setRandSeed(unsigned int);

  //! Allows the Coordinator class to add Global events to the
  //! Simulation.
  //!
  //! This is only allowed before the simulation is initialised.
  void addGlobal(Global*);

  //! Allows the Coordinator class to add System events to the
  //! Simulation.
  //!
  //! This is only allowed before the simulation is initialised.
  void addSystem(System*);

  //! Sets the ID of the simulation.
  //!
  //! This is used when running multiple simulations concurrently, and
  //! a ID is used to distinguish between them.
  //! \sa EReplicaExchangeSimulation
  //! \sa getSimID
  void setSimID(const size_t& n) { simID = n; }

  //! Returns the ID of the simulation.
  //! \sa setSimID
  size_t getSimID() { return simID; }

  //! Allows a Engine or the Coordinator class to access a named
  //! System event.
  //!
  //! This function is used by Engine classes like the
  //! EReplicaExchangeSimulation class to fetch System events like the
  //! Thermostat (e.g. CSysGhost).
  System* getSystem(std::string);

  //! Get the current time of the Simulation.
  //! \sa getnColl
  long double getSysTime();

  //! Get the Ensemble of the Simulation.
  inline const boost::scoped_ptr<dynamo::Ensemble>& getEnsemble() const 
  { return ensemble; }

  //! Get the Ensemble of the Simulation.
  inline boost::scoped_ptr<dynamo::Ensemble>& getEnsemble() 
  { return ensemble; }
  
  //! Get the std::ostringstream storing the Simulation history.
  //!
  //! This is used to add more lines to the Simulation's history,
  //! adding notes on what was done to create or run the Simulation.
  std::ostringstream &getHistory()
    { return ssHistory; }

  //! Get the number of Events executed.
  //! \sa getSysTime
  inline const unsigned long long &getnColl() const 
  { return eventCount; }

  //! Get the state of the Simulation.
  //! \sa ESimulationStatus
  inline const ESimulationStatus& getStatus() const
  { return status; }
  
  //! Allows a Engine or the Coordinator class to add an OutputPlugin
  //! to the Simulation.
  //! \param pluginDescriptor A string identifying a type of
  //! OutputPlugin to load and the values of any options to be
  //! parsed. E.g., "Plugin:OptA=1,OptB=2"
  void addOutputPlugin(std::string pluginDescriptor);
  
  //! Sets the frequency of the CSTicker event.
  void setTickerPeriod(double);

  //! Scales the frequency of the CSTicker event by the passed factor.
  void scaleTickerPeriod(double);
  
  /*! \brief An expensive sanity check for the system.
   *
   * Ensures that the system configuration is valid, and no overlaps
   * exist.
   */
  void checkSystem();
 private:
};

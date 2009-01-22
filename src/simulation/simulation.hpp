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
/*! \file simulation.hpp
 * \brief Contains the definition of the CSimulation class.
 */
#ifndef SIMULATION_H
#define SIMULATION_H

#include "../base/is_base.hpp"
#include "../base/is_simdata.hpp"
#include <boost/scoped_array.hpp>

class CDynamics;
class COutputPlugin;
class CStreamTask;

/*! \brief A single simulation holding particles, dynamics and output
 * plugins.
 *
 * This class is the typical realisation of a simulation program. It
 * can pretty much perform a standard simulation without any other
 * supporting class structure like the CEngine and CCoordinator. This
 * class handles the interface to the simulation and also stores the
 * CSimulation data by deriving from the SimData class.
 *
 *
 */
class CSimulation: public DYNAMO::Base_Class, public DYNAMO::SimData
{
 public:
  /*! \brief Just initialises the Base_Class
   */
  CSimulation();

  /*! \brief Initialise the entire CSimulation and the SimData struct.
   *
   * Most classes will have an initialisation function and its up to
   * this function to call them all and in the right order.
   */
  void initialise();

  inline void runSilentSimulation() { runSimulation(true); }
  //Gives a function signature to use
  inline void runSimulation() { runSimulation(false); }
  void runSimulation(bool silentMode);
  
  void loadXMLfile(const char *);
  
  void writeXMLfile(const char *);

  void loadPlugins(std::string);

  void outputData(const char* filename = "output.xml.bz2");
  
  void configLoaded();

  void simShutdown();

  void setTrajectoryLength(unsigned long long);

  void setnPrint(unsigned long long);

  void setnImage(unsigned long long);

  void setRandSeed(unsigned int);

  void addGlobal(CGlobal*);

  void addSystem(CSystem*);

  void initPlugins();

  void setSimID(const size_t& n) { simID = n; }

  CSystem* getSystem(std::string);

  long double getSysTime();

  inline const boost::scoped_ptr<DYNAMO::CEnsemble>& getEnsemble() const 
  { return Ensemble; }

  inline boost::scoped_ptr<DYNAMO::CEnsemble>& getEnsemble() 
  { return Ensemble; }
  
  std::ostringstream &getHistory()
    { return ssHistory; }

  inline const unsigned long long &getnColl() const 
  { return lNColl; }

  inline const ESimulationStatus& getStatus() const
  { return status; }

  inline void setBinaryXML(const bool& v) { binaryXML = v; }

  void addOutputPlugin(std::string);
  
  void setTickerPeriod(Iflt);

  void scaleTickerPeriod(Iflt);
  
 private:

  /*void executeEvent();

  void executeIntEvent();

  void executeGlobEvent();

  void executeLocalEvent();

  void executeSysEvent();*/
};

#endif

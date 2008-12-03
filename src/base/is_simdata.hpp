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

#ifndef IS_SimBase_H
#define IS_SimBase_H

#include "../dynamics/dynamics.hpp"
#include "../datatypes/vector.hpp"
#include "../simulation/particle.hpp"
#include "../outputplugins/outputplugin.hpp"
#include <vector>
#include <list>

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include "is_ensemble.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include "../extcode/include/boost/random/01_normal_distribution.hpp"

class CScheduler;
class CParticle;
class COutputPlugin;

template <class T>
class smrtPlugPtr;

//! \brief Holds the different phases of the simulation initialisation
typedef enum 
  {
    START         = 0, /*!< The first phase of the simulation. */
    CONFIG_LOADED = 1, /*!< After the configuration has been loaded. */
    INITIALISED   = 2, /*!< Once the classes have been initialised and
                          the simulation is ready to begin. */
    PRODUCTION    = 3, /*!< The simulation has already begun. */
    ERROR         = 4  /*!< The simulation has failed. */
  } ESimulationStatus;

namespace DYNAMO
{  
  typedef boost::mt19937 baseRNG;
  
  /*! \brief Fundamental collection of the CSimulation data.
   *
   * This struct contains all the data belonging to a single
   * CSimulation. It has been abstracted away from the CSimulation
   * class so that every class can contain a pointer to this datatype
   * without knowing the CSimulation class and causing a circular
   * reference/dependency.
   *
   * A pointer to this class has been incorporated in the base classes
   * SimBase and SimBase_Const which also provide some general
   * std::cout formatting.
   */
  struct SimData
  {
    /*! \brief Significant default value initialisation.
     */
    SimData();
    
    /*! \brief Handles deleting the CScheduler pointer
     *
     * \bug Make the Scheduler handled by a smrtPlugPtr, and make the
     * CScheduler class copyable.
     */
    ~SimData();
    
    /*! \brief Adds an output plugin of the type T.
     *
     * Throws an exception if there's already the same plugin loaded.
     */
    template<class T>
    const T* getOutputPlugin() const
    {
      BOOST_FOREACH(const smrtPlugPtr<COutputPlugin>& plugin, outputPlugins)
	if (dynamic_cast<const T*>(plugin.get_ptr()) != NULL)
	  return dynamic_cast<const T*>(plugin.get_ptr());
      
      D_throw() << "The output plugin " << (typeid(T).name()) << " is required, please add it";
    }   

    /*! \brief Finds a plugin of the given type using RTTI.
     *
     * Shouldn't use this frequently as its expensive, throws an
     * exception if it can't find the plugin.
     */
    template<class T>
    T* getOutputPlugin()
    {
      BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& plugin, outputPlugins)
	if (dynamic_cast<T*>(plugin.get_ptr()) != NULL)
	  return dynamic_cast<T*>(plugin.get_ptr());

      D_throw() << "The output plugin " << (typeid(T).name()) << " is required, please add it";
    }    

    /*! \brief The CEnsemble of the CSimulation. */
    boost::scoped_ptr<CEnsemble> Ensemble;

    /*! \brief The current system time of the simulation. 
     * 
     * This class is long double to reduce roundoff error as this gets
     * very large compared to an events delta t.
     */
    long double dSysTime;
    
    /*! \brief Number of events executed.*/
    unsigned long long lNColl;

    /*! \brief Maximum number of events to execute.*/
    unsigned long long lMaxNColl;

    /*! \brief How many events between periodic output/sampling*/
    unsigned long long lNPrint;
    
    /*! \brief Speeds the CSimulation loop by being the next periodic
        output collision number.*/
    unsigned long long lPrintLimiter;

    /*! \brief Number of CParticle's in the system. */
    unsigned long lN;
    
    /*! \brief The CParticle's of the system. */
    std::vector<CParticle> vParticleList;  
    
    /*! \brief A log of the previous simulation history. */
    std::ostringstream ssHistory;

    /*! \brief A ptr to the CScheduler of the system. */
    CScheduler *ptrScheduler;
    
    /*! \brief The CDynamics of the system. */
    CDynamics Dynamics;
    
    /*! \brief A vector of the ratio's of the simulation box/images sides.
     *
     * At least one ratio must be 1 as this is assumed when using the
     * ratio. i.e. it is normalised.
     */
    CVector<> aspectRatio;

    /*! \brief The random number generator of the system.
     */
    mutable baseRNG ranGenerator;

    mutable boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution_01<Iflt> > normal_sampler;

    mutable boost::uniform_01<DYNAMO::baseRNG, Iflt> uniform_sampler;  

    /*! \brief The collection of COutputPlugin's operating on this system.
     */
    std::vector<smrtPlugPtr<COutputPlugin> > outputPlugins; 

    /*! \brief The mean free time of the previous simulation run
     *
     * This is zero in the case that there is no previous simulation data.
     */
    Iflt lastRunMFT;

    /*! \brief This is just the ID number of the CSimulation in its container.
     *
     * This is used in the CEReplexer engine.
     */
    size_t simID;
    
    /*! \brief The current phase of the CSimulation.
     */
    ESimulationStatus status;
  };

}

#endif

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
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include "is_ensemble.hpp"

class CScheduler;
class CParticle;
class COutputPlugin;

template <class T>
class smrtPlugPtr;

typedef enum 
  {
    START         = 0,
    CONFIG_LOADED = 1,
    INITIALISED   = 2,
    PRODUCTION    = 3,
    ERROR         = 4
  } ESimulationStatus;

namespace DYNAMO
{  
  typedef boost::mt19937 baseRNG;
  
  struct SimData
  {     
    //Functions
    SimData();
    ~SimData();

    template<class T>
    const T* getOutputPlugin() const
    {
      BOOST_FOREACH(const smrtPlugPtr<COutputPlugin>& plugin, outputPlugins)
	if (dynamic_cast<const T*>(plugin.get_ptr()) != NULL)
	  return dynamic_cast<const T*>(plugin.get_ptr());

      I_throw() << "The output plugin " << (typeid(T).name()) << " is required, please add it";
    }    

    template<class T>
    T* getOutputPlugin()
    {
      BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& plugin, outputPlugins)
	if (dynamic_cast<T*>(plugin.get_ptr()) != NULL)
	  return dynamic_cast<T*>(plugin.get_ptr());

      I_throw() << "The output plugin " << (typeid(T).name()) << " is required, please add it";
    }    

    //Data structures
    boost::scoped_ptr<CEnsemble> Ensemble;
    long double dSysTime;
    unsigned long long lNColl, lMaxNColl, lNPrint, lPrintLimiter;
    unsigned long lN;
    size_t lReverseEvents;
    std::vector<CParticle> vParticleList;  
    std::ostringstream ssHistory;
    CScheduler *ptrScheduler;
    CDynamics Dynamics;
    CVector<> aspectRatio;
    mutable baseRNG ranGenerator;
    std::vector<smrtPlugPtr<COutputPlugin> > outputPlugins; 
    Iflt lastRunMFT; //The mean free time of the last simulation run
    unsigned int simID;
    ESimulationStatus status;
  };

}

#endif

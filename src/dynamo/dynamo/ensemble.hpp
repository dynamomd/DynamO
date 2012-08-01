/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#pragma once
#include <dynamo/base.hpp>
#include <tr1/array>

namespace magnet { namespace xml { class Node; class XmlStream; } }

namespace dynamo {
  class System;
  class Simulation;

  /*! \brief This class specifies the simulation ensemble that the
    simulation is being performed in.
   
    Sometimes, it is required to check the ensemble of the
    simulation. i.e. in replica exchange, we need to know if we're in
    the NVT ensemble, as this is the only time we can use the
    boltzmann relationship to calculate the exchange
    probability. Also, some plugins (thermal conductivity) are only
    valid in the NVE ensemble.

    The Ensemble for the simulated system is detected by the \ref
    Ensemble::loadEnsemble() function. Only the NVT and NVE Ensembles
    are defined, as these are the special cases we currently need to
    distinguish. More specialisations will be added as needed.
   */
  class Ensemble : public SimBase_const
  {
  public:
    Ensemble(const Simulation* const& SD, const char *aName = "Ensemble"):
      SimBase_const(SD, aName) {}
    
    virtual ~Ensemble() {}

    /*! \brief Used to determine what ensemble is correct for the
        current simulation.
      
      \param Sim A reference to the simulation data to detect the
      Ensemble for.
    */
    static shared_ptr<Ensemble> loadEnsemble(const dynamo::Simulation& Sim);

    /*! \brief Called to generate and store the Ensemble variables.
     */
    virtual void initialise() { dout << "Undefined Ensemble type."; }

    /*! \brief Returns an array containing the control values of the Ensemble
     (e.g., NVE) in the units of the output.

     \sa getEnsembleVals
    */
    virtual std::tr1::array<double,3> getReducedEnsembleVals() const { M_throw() << "Undefined Ensemble"; }
    
    /*! \brief Swaps the underlying ensemble control values.
      
      \sa EReplicaExchangeSimulation
    */
    virtual void swap(Ensemble& rhs) { std::swap(EnsembleVals, rhs.EnsembleVals); }

    /*! \brief Calculates the probability of carrying out a replica exchange
      move between this Ensemble and another.
    */
    virtual double exchangeProbability(const Ensemble&) const { M_throw() << "Undefined in this Ensemble"; }
    
    /*! Returns an array containing the ensemble values in simulation units.
    
      \sa getReducedEnsembleVals
    */
    virtual const std::tr1::array<double,3>& getEnsembleVals() const { M_throw() << "Undefined Ensemble"; }

  protected:
    std::tr1::array<double,3> EnsembleVals;
  };

  /*! \brief An Ensemble where N (no. of particles), V (simulation volume),
    and E (total energy) are held constant.
  */
  class EnsembleNVE           : public Ensemble 
  {
  public:
    EnsembleNVE(const dynamo::Simulation* SD): 
      Ensemble(SD, "EnsembleNVE") {}

    virtual void initialise();

    virtual std::tr1::array<double,3> getReducedEnsembleVals() const;
    
    virtual const std::tr1::array<double,3>& getEnsembleVals() const { return EnsembleVals; }
  };

  /*! \brief An Ensemble where N (no. of particles), V (simulation
   volume), and T (temperature) are held constant.
   
   This class also stores a pointer to the thermostat used to hold the
   temperature constant
  */
  class EnsembleNVT           : public Ensemble 
  {
  public:
    EnsembleNVT(const dynamo::Simulation* SD): 
      Ensemble(SD, "EnsembleNVT") {}

    virtual void initialise();

    virtual std::tr1::array<double,3> getReducedEnsembleVals() const;

    virtual double exchangeProbability(const Ensemble&) const;

    virtual const std::tr1::array<double,3>& getEnsembleVals() const { return EnsembleVals; }

  protected:
    shared_ptr<System> thermostat;
  };
}

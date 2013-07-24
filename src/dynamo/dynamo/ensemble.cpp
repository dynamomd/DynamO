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

#include <dynamo/ensemble.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/BC/PBC.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/dynamics/multicanonical.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  shared_ptr<Ensemble>
  Ensemble::loadEnsemble(const dynamo::Simulation& Sim)
  {
    bool hasThermostat = false; 

    try {
      std::shared_ptr<SysAndersen> thermoptr = std::dynamic_pointer_cast<SysAndersen>(Sim.systems["Thermostat"]);
      hasThermostat = bool(thermoptr);
    } catch (std::exception & err) {}
    
    //bool periodic = std::dynamic_pointer_cast<BCPeriodic>(Sim.BCs);
    
    if (hasThermostat)
      return shared_ptr<Ensemble>(new EnsembleNVT(&Sim));
    
    return shared_ptr<Ensemble>(new EnsembleNVE(&Sim));
    
    //Return the default unknown ensemble
    return shared_ptr<Ensemble>(new Ensemble(&Sim));
  }

  void
  EnsembleNVE::initialise()
  {
    EnsembleVals[0] = Sim->particles.size();
    EnsembleVals[1] = Sim->primaryCellSize[0] * Sim->primaryCellSize[1] * Sim->primaryCellSize[2];
    EnsembleVals[2] = Sim->calcInternalEnergy() + Sim->dynamics->getSystemKineticEnergy();

    dout << "NVE Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nV=" << EnsembleVals[1] / Sim->units.unitVolume()
	     << "\nE=" << EnsembleVals[2] / Sim->units.unitEnergy() << std::endl;
  }

  std::array<double,3> 
  EnsembleNVE::getReducedEnsembleVals() const
  {
    std::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->units.unitVolume();
    retval[2] = EnsembleVals[2] / Sim->units.unitEnergy();

    return retval;
  }

  void
  EnsembleNVT::initialise()
  {
    EnsembleVals[0] = Sim->particles.size();
    EnsembleVals[1] = Sim->units.unitVolume();

    try {
      thermostat = Sim->systems["Thermostat"];
    } catch (std::exception & err)
      {
	M_throw() << "Could not find the Thermostat in NVT system\n"
		  << err.what();
      }
    
    //Only one kind of thermostat so far!
    if (!std::dynamic_pointer_cast<SysAndersen>(thermostat))
      {
	M_throw() << "Could not upcast thermostat to Andersens";
      }    
    
    EnsembleVals[2] = std::dynamic_pointer_cast<SysAndersen>(thermostat)->getTemperature();
    
    dout << "NVT Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nV=" << EnsembleVals[1] / Sim->units.unitVolume()
	     << "\nT=" << EnsembleVals[2] / Sim->units.unitEnergy() << std::endl;
  }

  std::array<double,3> 
  EnsembleNVT::getReducedEnsembleVals() const
  {
    std::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->units.unitVolume();
    retval[2] = EnsembleVals[2] / Sim->units.unitEnergy();

    return retval;
  }

  double 
  EnsembleNVT::exchangeProbability(const Ensemble& oE) const
  {
#ifdef DYNAMO_DEBUG
    if (dynamic_cast<const EnsembleNVT*>(&oE) == NULL)
      M_throw() << "The ensembles types differ";
#endif

    //Must use static cast to allow access to protected members
    
    const EnsembleNVT& ensemble2(static_cast<const EnsembleNVT&>(oE));

    double beta1 = 1 / EnsembleVals[2];
    double E1 = Sim->getOutputPlugin<OPMisc>()->getConfigurationalU();
    double beta2 = 1 / ensemble2.getEnsembleVals()[2];
    double E2 = ensemble2.Sim->getOutputPlugin<OPMisc>()->getConfigurationalU();
    
    //This is -\Delta in the Sugita_Okamoto paper
    double factor = (E1 - E2) * (beta1 - beta2);
    
    if (std::dynamic_pointer_cast<DynNewtonianMC>(Sim->dynamics))
      {
	factor += static_cast<const DynNewtonianMC&>(*Sim->dynamics).W(E1);
	factor -= static_cast<const DynNewtonianMC&>(*Sim->dynamics).W(E2);
      }

    if (std::dynamic_pointer_cast<DynNewtonianMC>(ensemble2.Sim->dynamics))
      {
	factor += static_cast<const DynNewtonianMC&>(*ensemble2.Sim->dynamics).W(E2);
	factor -= static_cast<const DynNewtonianMC&>(*ensemble2.Sim->dynamics).W(E1);
      }

    return std::exp(factor);
  }
}

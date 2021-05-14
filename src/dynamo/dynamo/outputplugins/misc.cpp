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

#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <magnet/memUsage.hpp>
#include <magnet/xmlwriter.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <ctime>

namespace dynamo {
  OPMisc::OPMisc(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"Misc",0),//ContactMap must be after this
    _dualEvents(0),
    _singleEvents(0),
    _virtualEvents(0),
    _reverseEvents(0),
    _lateInitComplete(false)
  {}

  void
  OPMisc::replicaExchange(OutputPlugin& misc2)
  {
    //We must swap anything that is associated with the sampling
    //(averages, sums) but keep anything that is related to the
    //configuration at this instant of time.

    OPMisc& op = static_cast<OPMisc&>(misc2);
    
    std::swap(_counters, op._counters);
    std::swap(_starttime, op._starttime);
    std::swap(_dualEvents, op._dualEvents);
    std::swap(_singleEvents, op._singleEvents);
    std::swap(_virtualEvents, op._virtualEvents);
    std::swap(_reverseEvents, op._reverseEvents);
    
    _KE.swapAverages(op._KE);
    _internalE.swapAverages(op._internalE);
    _sysMomentum.swapAverages(op._sysMomentum);
    _kineticP.swapAverages(op._kineticP);

    std::swap(collisionalP, op.collisionalP);

    _thermalConductivity.clear();
    _viscosity.clear();
    _bulkVisc.clear();
    _crossVisc.clear();

    for (auto& correlator : _thermalDiffusion)
      correlator.clear();
    for (auto& correlator : _mutualDiffusion)
      correlator.clear();

    //These remain unchanged
    //_internalEnergy;
    //_speciesMasses;
    //_speciesMomenta;
    //_systemMass;
  }

  void
  OPMisc::temperatureRescale(const double& scale)
  { 
    _KE  = _KE.current() * scale;
  }

  double 
  OPMisc::getMeankT() const
  {
    return 2.0 * _KE.mean() / Sim->dynamics->getParticleDOF();
  }

  double 
  OPMisc::getMeanSqrkT() const
  {
    return 4.0 * _KE.meanSqr() / std::pow(Sim->dynamics->getParticleDOF(), 2);
  }

  double 
  OPMisc::getCurrentkT() const
  {
    return 2.0 * _KE.current() / Sim->dynamics->getParticleDOF();
  }

  double 
  OPMisc::getMeanUConfigurational() const
  { 
    return _internalE.mean(); 
  }

  double 
  OPMisc::getMeanSqrUConfigurational() const
  { return _internalE.meanSqr(); }

  void outputCorrelator(magnet::xml::XmlStream& XML, const bool i1, const bool i2, const double inv_units, const double time_units, magnet::math::LogarithmicTimeCorrelator<double>& corr) {
    std::string t = "CC";
    if (i1) t[0] = 'I';
    if (i2) t[1] = 'I';

    using namespace magnet::xml;
    
    XML	<< tag("Component") << attr("type") << t
	<< chardata();
    {
      std::vector<magnet::math::LogarithmicTimeCorrelator<double>::Data>
	data = corr.getAveragedCorrelator(i1, i2);
            
      XML << "0 0 0\n";
      for (size_t i(0); i < data.size(); ++i) {
	XML << data[i].time / time_units << " "
	    << data[i].sample_count << " "
	    << data[i].value * inv_units << " "
	    << "\n";
      }
    }

    XML << endtag("Component");
  }

  void outputCorrelator(magnet::xml::XmlStream& XML, const bool i1, const bool i2, const double inv_units, const double time_units, magnet::math::LogarithmicTimeCorrelator<Vector>& corr) {
    std::string t = "CC";
    if (i1) t[0] = 'I';
    if (i2) t[1] = 'I';

    using namespace magnet::xml;
    
    XML	<< tag("Component") << attr("type") << t
	<< chardata();
    {
      std::vector<magnet::math::LogarithmicTimeCorrelator<Vector>::Data>
	data = corr.getAveragedCorrelator(i1, i2);
            
      XML << "0 0 0 0 0\n";
      for (size_t i(0); i < data.size(); ++i) {
	XML << data[i].time / time_units << " "
	    << data[i].sample_count << " ";
	for (size_t j(0); j < 3; ++j)
	  XML << data[i].value[j] * inv_units << " ";
	XML << "\n";
      }
    }

    XML << endtag("Component");
  }
  
  void outputCorrelator(magnet::xml::XmlStream& XML, const bool i1, const bool i2, const double inv_units, const double time_units, magnet::math::LogarithmicTimeCorrelator<Matrix>& corr) {
    std::string t = "CC";
    if (i1) t[0] = 'I';
    if (i2) t[1] = 'I';

    using namespace magnet::xml;
    
    XML	<< tag("Component") << attr("type") << t
	<< chardata();
    
    {
      std::vector<magnet::math::LogarithmicTimeCorrelator<Matrix>::Data>
	data = corr.getAveragedCorrelator(i1, i2);
      

      XML << "0 0 0 0 0 0 0 0 0 0 0\n";
      for (size_t i(0); i < data.size(); ++i)
	{
	  XML << data[i].time / time_units << " "
	      << data[i].sample_count << " ";
	  
	  for (size_t j(0); j < 3; ++j)
	    for (size_t k(0); k < 3; ++k)
	      XML << data[i].value(j, k) * inv_units << " ";
	  XML << "\n";
	}
    }
    
    XML << endtag("Component");    
  }

  template<class T>
  void outputCorrelator(magnet::xml::XmlStream& XML, const double inv_units, const double time_units, magnet::math::LogarithmicTimeCorrelator<T>& corr) {
    outputCorrelator(XML, false, false, inv_units, time_units, corr);
    outputCorrelator(XML, false, true, inv_units, time_units, corr);
    outputCorrelator(XML, true, false, inv_units, time_units, corr);
    outputCorrelator(XML, true, true, inv_units, time_units, corr);    
  }

  void
  OPMisc::initialise()
  {
    _KE.init(Sim->dynamics->getSystemKineticEnergy());
    _internalE.init(Sim->calcInternalEnergy());

    dout << "Particle Count " << Sim->N()
	 << "\nSim Unit Length " << Sim->units.unitLength()
	 << "\nSim Unit Time " << Sim->units.unitTime()
	 << "\nDensity " << Sim->getNumberDensity() * Sim->units.unitVolume()
	 << "\nPacking Fraction " << Sim->getPackingFraction()
	 << "\nTemperature " << getCurrentkT() / Sim->units.unitEnergy();
    
    if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs)) {
      dout << " (Assuming linear shear profile)" <<std::endl;
      derr << "\nTemperature output disabled, as shearing/LEBC can give non-linear velocity profiles." << std::endl;
    }
    
    dout << "\nNo. of Species " << Sim->species.size()
	 << "\nSimulation box length "
	 << (Sim->primaryCellSize / Sim->units.unitLength()).toString()
	 << std::endl;

    Matrix kineticP;
    Vector thermalConductivityFS({0, 0, 0});
    _speciesMomenta.clear();
    _speciesMomenta.resize(Sim->species.size());
    _speciesMasses.clear();
    _speciesMasses.resize(Sim->species.size());

    _internalEnergy.clear();
    _internalEnergy.resize(Sim->N(), 0);

    for (const auto& p1 : Sim->particles)
      {
	std::unique_ptr<IDRange> ids(Sim->ptrScheduler->getParticleNeighbours(p1));
	for (size_t ID2 : *ids)
	  if (ID2 != p1.getID())
	    _internalEnergy[p1.getID()] += 0.5 * Sim->getInteraction(p1, Sim->particles[ID2])->getInternalEnergy(p1, Sim->particles[ID2]);
      }

    for (const Particle& part : Sim->particles)
      {
	const Species& sp = *(Sim->species(part));
	const double mass = sp.getMass(part.getID());
	if (std::isinf(mass)) continue;
	kineticP += mass * Dyadic(part.getVelocity(), part.getVelocity());
	_speciesMasses[sp.getID()] += mass;
	_speciesMomenta[sp.getID()] += mass * part.getVelocity();
	thermalConductivityFS += part.getVelocity() * (sp.getParticleKineticEnergy(part) + _internalEnergy[part.getID()]);
      }

    Vector sysMomentum({0, 0, 0});
    _systemMass = 0;
    for (size_t i(0); i < Sim->species.size(); ++i)
      {
	sysMomentum += _speciesMomenta[i];
	_systemMass += _speciesMasses[i];
      }

    _kineticP.init(kineticP);
    _sysMomentum.init(sysMomentum);

    //Set up the correlators
    double correlator_dt = Sim->lastRunMFT / 8;
    if (correlator_dt == 0.0)
      correlator_dt = 1.0 / sqrt(getCurrentkT());

    _thermalConductivity.resize(correlator_dt, 10, 2, false);
    _thermalConductivity.setFreeStreamValue(thermalConductivityFS);

    _viscosity.resize(correlator_dt, 10, 2, true);
    _viscosity.setFreeStreamValue(kineticP);
    
    _bulkVisc.resize(correlator_dt, 10, 2, true);
    double isoViscFS(0);
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      isoViscFS += kineticP(iDim, iDim);
    _bulkVisc.setFreeStreamValue(isoViscFS / 3);
    
    _crossVisc.resize(correlator_dt, 10, 2, true);
    Vector crossViscFS1({0, 0, 0});
    Vector crossViscFS2({0, 0, 0});
    for (size_t iDim(0); iDim < NDIM; ++iDim) {
      crossViscFS1[iDim] = kineticP(iDim, iDim);
      crossViscFS2[iDim] = kineticP((iDim+1) % NDIM, (iDim+1) % NDIM);
    }
    _crossVisc.setFreeStreamValue(crossViscFS1, crossViscFS2);

    _thermalDiffusion.resize(Sim->species.size());
    _mutualDiffusion.resize(Sim->species.size() * Sim->species.size());
    for (size_t spid1(0); spid1 < Sim->species.size(); ++spid1)
      {
	_thermalDiffusion[spid1].resize(correlator_dt, 10, 2, false);
	_thermalDiffusion[spid1].setFreeStreamValue(thermalConductivityFS, _speciesMomenta[spid1] - (_speciesMasses[spid1] / _systemMass) * _sysMomentum.current());
	
	for (size_t spid2(spid1); spid2 < Sim->species.size(); ++spid2)
	  {
	    _mutualDiffusion[spid1 * Sim->species.size() + spid2].resize(correlator_dt, 10, 2, false);
	    _mutualDiffusion[spid1 * Sim->species.size() + spid2].setFreeStreamValue
	      (_speciesMomenta[spid1] - (_speciesMasses[spid1] / _systemMass) * _sysMomentum.current(),
	       _speciesMomenta[spid2] - (_speciesMasses[spid2] / _systemMass) * _sysMomentum.current());
	  }
      }

    dout << "Total momentum < ";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << _sysMomentum.current()[iDim] / Sim->units.unitMomentum() << " ";
    dout << ">" << std::endl;

    _starttime = std::chrono::system_clock::now();
  }

  void
  OPMisc::eventUpdate(const Event& eevent, const NEventData& NDat)
  {
    //if ((!_lateInitComplete) && (Sim->eventCount > 100 * Sim->N())) {
    //  _lateInitComplete = true;
    //  const auto avgPk = _kineticP.mean();
    //  dout << "## Reinit using Pk(0,0)=" << avgPk(0,0) << std::endl; 
    //  const auto avgPc = collisionalP / (Sim->systemTime - eevent._dt);
    //  dout << "## Reinit using Pc(0,0)=" << avgPc(0,0) << std::endl; 
    //  _viscosity.resize(_viscosity.getSampleTime(), 10, 2, avgPk, avgPk, avgPc, avgPc);
    //}
    
    stream(eevent._dt);
    auto type = eevent._type;

    if ((NDat.L1partChanges.size() == 1) && (NDat.L2partChanges.size() == 0))
      type = NDat.L1partChanges[0].getType();
    if ((NDat.L1partChanges.size() == 0) && (NDat.L2partChanges.size() == 1))
      type = NDat.L2partChanges[0].getType();
    
    CounterData& counterdata = _counters[CounterKey(getClassKey(eevent), type)];
    counterdata.count += NDat.L1partChanges.size() + NDat.L2partChanges.size();

    Vector thermalDel({0,0,0});
    for (const ParticleEventData& PDat : NDat.L1partChanges)
      {

	const Particle& part = Sim->particles[PDat.getParticleID()];
	const Species& species = *Sim->species(part);
	const double mass = species.getMass(part.getID());
	const double deltaKE = species.getParticleKineticEnergy(part) - PDat.getOldKE();
	
	_KE += deltaKE;
	_internalE += PDat.getDeltaU();
	//This must be updated before p1E is calculated
	_internalEnergy[PDat.getParticleID()] += PDat.getDeltaU();
	const double p1E = species.getParticleKineticEnergy(part) + _internalEnergy[PDat.getParticleID()];
	const double p1deltaE = deltaKE + PDat.getDeltaU();
	Vector delP1 = mass * (part.getVelocity() - PDat.getOldVel());

	counterdata.netimpulse += mass * (part.getVelocity() -  PDat.getOldVel());
	counterdata.netKEchange += deltaKE;
	counterdata.netUchange += PDat.getDeltaU();

        _singleEvents += (PDat.getType() != VIRTUAL);
	_virtualEvents += (PDat.getType() == VIRTUAL);
	
	_kineticP += mass * (Dyadic(part.getVelocity(), part.getVelocity()) - Dyadic(PDat.getOldVel(), PDat.getOldVel()));
	_sysMomentum += delP1;
	_speciesMomenta[species.getID()] += delP1;
	thermalDel += part.getVelocity() * p1E - PDat.getOldVel() * (p1E - p1deltaE);
      }

    for (const PairEventData& PDat : NDat.L2partChanges)
      {
	const Particle& part1 = Sim->particles[PDat.particle1_.getParticleID()];
	const Particle& part2 = Sim->particles[PDat.particle2_.getParticleID()];
	const Species& sp1 = *Sim->species[PDat.particle1_.getSpeciesID()];
	const Species& sp2 = *Sim->species[PDat.particle2_.getSpeciesID()];

        const double p1KE = sp1.getParticleKineticEnergy(part1);
	const double p1E = p1KE + _internalEnergy[part1];
        
        const double p2KE = sp2.getParticleKineticEnergy(part2);
	const double p2E = p2KE + _internalEnergy[part2];

	const double deltaKE1 = p1KE - PDat.particle1_.getOldKE();
	const double deltaKE2 = p2KE - PDat.particle2_.getOldKE();
	const double p1deltaE = deltaKE1 + PDat.particle1_.getDeltaU();
        const double p2deltaE = deltaKE2 + PDat.particle2_.getDeltaU();
        
	const double mass1 = sp1.getMass(part1.getID());
	const double mass2 = sp2.getMass(part2.getID());
	const Vector delP = mass1 * (part1.getVelocity() - PDat.particle1_.getOldVel());


	const double deltaKE = deltaKE1 + deltaKE2;
	_KE += deltaKE;
	counterdata.netKEchange += deltaKE;
	const double deltaU = PDat.particle1_.getDeltaU() + PDat.particle2_.getDeltaU();
	_internalE += deltaU;
	counterdata.netUchange += deltaU;

	_internalEnergy[PDat.particle1_.getParticleID()] += PDat.particle1_.getDeltaU();
	_internalEnergy[PDat.particle2_.getParticleID()] += PDat.particle2_.getDeltaU();
	_dualEvents += (PDat.getType() != VIRTUAL);
	_virtualEvents += (PDat.getType() == VIRTUAL);

	collisionalP += magnet::math::Dyadic(PDat.rij, delP);

	_kineticP
	  += mass1 * (Dyadic(part1.getVelocity(), part1.getVelocity())
		      - Dyadic(PDat.particle1_.getOldVel(), PDat.particle1_.getOldVel()))
	  + mass2 * (Dyadic(part2.getVelocity(), part2.getVelocity())
		     - Dyadic(PDat.particle2_.getOldVel(), PDat.particle2_.getOldVel()));

	const auto visc_imp = magnet::math::Dyadic(PDat.rij, delP);
	_viscosity.addImpulse(visc_imp);


	double isoVisc_imp(0);
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  isoVisc_imp += visc_imp(iDim, iDim);
	_bulkVisc.addImpulse(isoVisc_imp / 3);
    
	Vector crossVisc_imp1({0, 0, 0});
	Vector crossVisc_imp2({0, 0, 0});
	for (size_t iDim(0); iDim < NDIM; ++iDim) {
	  crossVisc_imp1[iDim] = visc_imp(iDim, iDim);
	  crossVisc_imp2[iDim] = visc_imp((iDim+1) % NDIM, (iDim+1) % NDIM);
	}
	_crossVisc.addImpulse(crossVisc_imp1, crossVisc_imp2);
	
	
	_speciesMomenta[sp1.getID()] += delP;
	_speciesMomenta[sp2.getID()] -= delP;

	const Vector thermalImpulse = PDat.rij * p1deltaE;

	_thermalConductivity.addImpulse(thermalImpulse);

	for (size_t spid1(0); spid1 < Sim->species.size(); ++spid1)
	  _thermalDiffusion[spid1].addImpulse(thermalImpulse, Vector{0,0,0});

	thermalDel += part1.getVelocity() * p1E + part2.getVelocity() * p2E
	  - PDat.particle1_.getOldVel() * (p1E - p1deltaE) - PDat.particle2_.getOldVel() * (p2E - p2deltaE);
      }

    _thermalConductivity.setFreeStreamValue
      (_thermalConductivity.getFreeStreamValue() + thermalDel);

    const auto kineticP = _kineticP.current();
    
    _viscosity.setFreeStreamValue(kineticP);

    double isoViscFS(0);
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      isoViscFS += kineticP(iDim, iDim);
    _bulkVisc.setFreeStreamValue(isoViscFS / 3);
    
    Vector crossViscFS1({0, 0, 0});
    Vector crossViscFS2({0, 0, 0});
    for (size_t iDim(0); iDim < NDIM; ++iDim) {
      crossViscFS1[iDim] = kineticP(iDim, iDim);
      crossViscFS2[iDim] = kineticP((iDim+1) % NDIM, (iDim+1) % NDIM);
    }
    _crossVisc.setFreeStreamValue(crossViscFS1, crossViscFS2);

    
    for (size_t spid1(0); spid1 < Sim->species.size(); ++spid1)
      {
	_thermalDiffusion[spid1]
	  .setFreeStreamValue(_thermalConductivity.getFreeStreamValue(),
			      _speciesMomenta[spid1] - _sysMomentum.current() * (_speciesMasses[spid1] / _systemMass));

	for (size_t spid2(spid1); spid2 < Sim->species.size(); ++spid2)
	  _mutualDiffusion[spid1 * Sim->species.size() + spid2].setFreeStreamValue
	    (_speciesMomenta[spid1] - (_speciesMasses[spid1] / _systemMass) * _sysMomentum.current(),
	     _speciesMomenta[spid2] - (_speciesMasses[spid2] / _systemMass) * _sysMomentum.current());
      }
  }

  void
  OPMisc::stream(double dt)
  {
    _reverseEvents += (dt < 0);
    _KE.stream(dt);
    _internalE.stream(dt);
    _kineticP.stream(dt);
    _sysMomentum.stream(dt);
    _thermalConductivity.freeStream(dt);
    _viscosity.freeStream(dt);
    _bulkVisc.freeStream(dt);
    _crossVisc.freeStream(dt);
    for (size_t spid1(0); spid1 < Sim->species.size(); ++spid1)
      {
	_thermalDiffusion[spid1].freeStream(dt);
	for (size_t spid2(spid1); spid2 < Sim->species.size(); ++spid2)
	  _mutualDiffusion[spid1 * Sim->species.size() + spid2].freeStream(dt);
      }
  }

  double
  OPMisc::getMFT() const
  {
    return Sim->systemTime * static_cast<double>(Sim->N())
      /(Sim->units.unitTime() * ((2.0 * static_cast<double>(_dualEvents)) + static_cast<double>(_singleEvents)));
  }

  double 
  OPMisc::getDuration() const {
    return std::chrono::duration<double>(std::chrono::system_clock::now() - _starttime).count();
  }
  
  double 
  OPMisc::getEventsPerSecond() const {
    return Sim->eventCount / getDuration();
  }

  double 
  OPMisc::getSimTimePerSecond() const {
    return Sim->systemTime / (getDuration() * Sim->units.unitTime());
  }
  
  Matrix 
  OPMisc::getPressureTensor() const
  { return ((collisionalP / Sim->systemTime) + _kineticP.mean()) / Sim->getSimVolume(); }

  void
  OPMisc::output(magnet::xml::XmlStream &XML)
  {
    using namespace magnet::xml;

    dout << "\nTotal events executed " << Sim->eventCount
	 << "\nSimulation end time  " << Sim->systemTime / Sim->units.unitTime()
	 << "\nAvg. events/s " << getEventsPerSecond()
	 << "\nSim time per second " << getSimTimePerSecond()
	 << std::endl;

    const double V = Sim->getSimVolume();
    const Matrix collP = collisionalP / (V * Sim->systemTime);
    const Matrix P = (_kineticP.mean() / V) + collP;

    XML << tag("Misc")

	<< tag("Timing")
	<< attr("RuntimeSeconds") <<  getDuration()
	<< attr("RuntimeHours") <<  getDuration() / 3600
	<< attr("EventsPerSec") << getEventsPerSecond()
	<< attr("SimTimePerSec") << getSimTimePerSecond()
	<< endtag("Timing")

	<< tag("Density")
	<< attr("val")
	<< Sim->getNumberDensity() * Sim->units.unitVolume()
	<< endtag("Density")

	<< tag("PackingFraction")
	<< attr("val") << Sim->getPackingFraction()
	<< endtag("PackingFraction")

	<< tag("SpeciesCount")
	<< attr("val") << Sim->species.size()
	<< endtag("SpeciesCount")

	<< tag("ParticleCount")
	<< attr("val") << Sim->N()
	<< endtag("ParticleCount")

	<< tag("SystemMomentum")
	<< tag("Current")
	<< attr("x") << _sysMomentum.current()[0] / Sim->units.unitMomentum()
	<< attr("y") << _sysMomentum.current()[1] / Sim->units.unitMomentum()
	<< attr("z") << _sysMomentum.current()[2] / Sim->units.unitMomentum()
	<< endtag("Current")
	<< tag("Average")
	<< attr("x") << _sysMomentum.mean()[0] / Sim->units.unitMomentum()
	<< attr("y") << _sysMomentum.mean()[1] / Sim->units.unitMomentum()
	<< attr("z") << _sysMomentum.mean()[2] / Sim->units.unitMomentum()
	<< endtag("Average")
	<< endtag("SystemMomentum");

    if (!std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      XML << tag("Temperature")
	  << attr("Mean") << getMeankT() / Sim->units.unitEnergy()
	  << attr("MeanSqr") << getMeanSqrkT() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	  << attr("Current") << getCurrentkT() / Sim->units.unitEnergy()
	  << attr("Min") << 2.0 * _KE.min() / (Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	  << attr("Max") << 2.0 * _KE.max() / (Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	  << endtag("Temperature");
    
    XML << tag("UConfigurational")
	<< attr("Mean") << getMeanUConfigurational() / Sim->units.unitEnergy()
	<< attr("MeanSqr") << getMeanSqrUConfigurational() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	<< attr("Current") << _internalE.current() / Sim->units.unitEnergy()
	<< attr("Min") << _internalE.min() / Sim->units.unitEnergy()
	<< attr("Max") << _internalE.max() / Sim->units.unitEnergy()
	<< endtag("UConfigurational");
    
    if (!std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      XML << tag("ResidualHeatCapacity")
	  << attr("Value") 
	  << (getMeanSqrUConfigurational() - getMeanUConfigurational() * getMeanUConfigurational())
	/ (getMeankT() * getMeankT())
	  << endtag("ResidualHeatCapacity");
    
    XML << tag("Pressure")
	<< attr("Avg") << P.tr() / (3.0 * Sim->units.unitPressure())
	<< tag("Tensor") << chardata()
      ;
    
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  XML << P(iDim, jDim) / Sim->units.unitPressure() << " ";
	XML << "\n";
      }
    
    XML << endtag("Tensor")
	<< tag("InteractionContribution") << chardata()
      ;
    
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  XML << collP(iDim, jDim) / Sim->units.unitPressure() << " ";
	XML << "\n";
      }
    
    XML << endtag("InteractionContribution")
	<< endtag("Pressure")

	<< tag("Duration")
	<< attr("Events") << Sim->eventCount
	<< attr("OneParticleEvents") << _singleEvents
	<< attr("TwoParticleEvents") << _dualEvents
	<< attr("VirtualEvents") << _virtualEvents
	<< attr("Time") << Sim->systemTime / Sim->units.unitTime()
	<< endtag("Duration")

	<< tag("EventCounters");
  
    typedef std::pair<CounterKey, CounterData> mappair;
    for (const mappair& mp1 : _counters)
      XML << tag("Entry")
	  << attr("Type") << getClass(mp1.first.first)
	  << attr("Name") << getName(mp1.first.first, Sim)
	  << attr("Event") << mp1.first.second
	  << attr("Count") << mp1.second.count
	  << tag("NetImpulse") 
	  << mp1.second.netimpulse / Sim->units.unitMomentum()
	  << endtag("NetImpulse")
	  << tag("NetKEChange")
	  << attr("Value") << mp1.second.netKEchange / Sim->units.unitEnergy()
	  << endtag("NetKEChange")
	  << tag("NetUChange")
	  << attr("Value") << mp1.second.netUchange / Sim->units.unitEnergy()
	  << endtag("NetUChange")
	  << endtag("Entry");
    
    XML << endtag("EventCounters")

	<< tag("PrimaryImageSimulationSize")
	<< Sim->primaryCellSize / Sim->units.unitLength()
	<< endtag("PrimaryImageSimulationSize")

	<< tag("totMeanFreeTime")
	<< attr("val")
	<< getMFT()
	<< endtag("totMeanFreeTime")

	<< tag("NegativeTimeEvents")
	<< attr("Count") << _reverseEvents
	<< endtag("NegativeTimeEvents")

	<< tag("Memusage")
	<< attr("MaxKiloBytes") << magnet::process_mem_usage()
	<< endtag("Memusage");

    if (!std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs)) {
      XML << tag("ThermalConductivity") << tag("Correlator");

      {
	const double inv_units = Sim->units.unitk()
	  / ( Sim->units.unitTime() * Sim->units.unitThermalCond() * 2.0 * getMeankT() * V);
	
	outputCorrelator(XML, inv_units, Sim->units.unitTime(), _thermalConductivity);
      }
      
      XML << endtag("Correlator") << endtag("ThermalConductivity")
	  << tag("Viscosity") << tag("Correlator");

      {
	const double inv_units = 1.0 / (Sim->units.unitTime() * Sim->units.unitViscosity() * 2.0 * getMeankT() * V);
	outputCorrelator(XML, inv_units, Sim->units.unitTime(), _viscosity);
      }
      
      XML << endtag("Correlator") << endtag("Viscosity")
	  << tag("BulkViscosity") << tag("Correlator");

      {
	const double inv_units = 1.0 / (Sim->units.unitTime() * Sim->units.unitViscosity() * 2.0 * getMeankT() * V);
	outputCorrelator(XML, inv_units, Sim->units.unitTime(), _bulkVisc);
      }
      
      XML << endtag("Correlator") << endtag("BulkViscosity")
	  << tag("CrossViscosity") << tag("Correlator");

      {
	const double inv_units = 1.0 / (Sim->units.unitTime() * Sim->units.unitViscosity() * 2.0 * getMeankT() * V);
	outputCorrelator(XML, inv_units, Sim->units.unitTime(), _crossVisc);
      }
      
      XML << endtag("Correlator") << endtag("CrossViscosity")
	  << tag("ThermalDiffusion");
      
      for (size_t i(0); i < Sim->species.size(); ++i)
	{
	  XML << tag("Correlator")
	      << attr("Species") << Sim->species[i]->getName();

	  const double inv_units = 1.0
	    / (Sim->units.unitTime() * Sim->units.unitThermalDiffusion() * 2.0 * getMeankT() * V);
	  
	  outputCorrelator(XML, inv_units, Sim->units.unitTime(), _thermalDiffusion[i]);
	  
	  XML << endtag("Correlator");
	}

      XML << endtag("ThermalDiffusion")
	  << tag("MutualDiffusion");

      for (size_t i(0); i < Sim->species.size(); ++i)
	for (size_t j(i); j < Sim->species.size(); ++j)
	  {
	    XML << tag("Correlator")
		<< attr("Species1") << Sim->species[i]->getName()
		<< attr("Species2") << Sim->species[j]->getName();

	    const double inv_units = 1.0
	      / (Sim->units.unitTime() * Sim->units.unitMutualDiffusion() * 2.0 * getMeankT() * V);
	    
	    outputCorrelator(XML, inv_units, Sim->units.unitTime(), _mutualDiffusion[i * Sim->species.size() + j]);

	    XML << endtag("Correlator");
	  }

      XML << endtag("MutualDiffusion");
    }
    XML << endtag("Misc");
  }

  void
  OPMisc::periodicOutput()
  {
    //Calculate the ETA of the simulation, and take care with overflows and the like
    double _earliest_end_time = std::numeric_limits<float>::infinity();
    for (const auto& sysPtr : Sim->systems)
      if (std::dynamic_pointer_cast<SystHalt>(sysPtr))
	_earliest_end_time = std::min(_earliest_end_time, sysPtr->getdt());

    double time_seconds_remaining = _earliest_end_time / (getSimTimePerSecond() * Sim->units.unitTime());
    size_t seconds_remaining = time_seconds_remaining;
    
    if (time_seconds_remaining > std::numeric_limits<size_t>::max())
      seconds_remaining = std::numeric_limits<size_t>::max();

    if (Sim->endEventCount != std::numeric_limits<size_t>::max())
      {
	double event_seconds_remaining = (Sim->endEventCount - Sim->eventCount) / getEventsPerSecond() + 0.5;
	
	if (event_seconds_remaining < std::numeric_limits<size_t>::max())
	  seconds_remaining = std::min(seconds_remaining, size_t(event_seconds_remaining));

      }
  
    if (seconds_remaining != std::numeric_limits<size_t>::max())
      {
	size_t ETA_days = seconds_remaining / (3600 * 24);
	size_t ETA_hours = (seconds_remaining / 3600) % 24;
	size_t ETA_mins = (seconds_remaining / 60) % 60;
	size_t ETA_secs = seconds_remaining % 60;

	//Here we try to give the friendliest ETA, as people can get
	//really tense about simulation timings. We only show two
	//places (hrs/min) (min/s) as the error in the estimate can be
	//substantial, particularly at the start.
	I_Pcout() << "ETA ";
	if (ETA_days)
	  //Show days, and round to nearest hour
	  I_Pcout() << ETA_days << "d " << std::lround(ETA_hours + float(ETA_mins) / 60) << "hr, ";
	else if (ETA_hours)
	  //Show hours, and round to nearest minute
	  I_Pcout() << ETA_hours << "hr " << std::lround(ETA_mins + float(ETA_secs) / 60) << "min, ";
	else if (ETA_mins > 5)
	  I_Pcout() << std::lround(ETA_mins + float(ETA_secs) / 60) << "min, ";
	else if (ETA_mins)
	  I_Pcout() << ETA_mins << "min " << ETA_secs << "s, ";
	else
	  I_Pcout() << ETA_secs << "s, ";
      
    }
    
    I_Pcout() << "Events " << (Sim->eventCount+1)/1000 << "k, t "
	      << Sim->systemTime/Sim->units.unitTime() 
	      << ", <MFT> " <<  getMFT();
    
    if (!std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      I_Pcout() << ", T " << getCurrentkT() / Sim->units.unitEnergy();
    
    I_Pcout() << ", U " << _internalE.current() / (Sim->units.unitEnergy() * Sim->N());
  }
}

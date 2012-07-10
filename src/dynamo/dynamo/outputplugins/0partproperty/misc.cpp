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

#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/memUsage.hpp>
#include <magnet/xmlwriter.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <boost/foreach.hpp>
#include <sys/time.h>
#include <ctime>

namespace dynamo {
  OPMisc::OPMisc(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"Misc",0),
    _dualEvents(0),
    _singleEvents(0),
    _virtualEvents(0),
    _reverseEvents(0)
  {}

  void
  OPMisc::changeSystem(OutputPlugin* misc2)
  {
    OPMisc& op = static_cast<OPMisc&>(*misc2);
    _KE.swapCurrentValues(op._KE);
    _internalE.swapCurrentValues(op._internalE);
    _kineticP.swapCurrentValues(op._kineticP);

    std::swap(Sim, op.Sim);
  }

  void
  OPMisc::temperatureRescale(const double& scale)
  { 
    _KE  = _KE.current() * scale;
  }

  double 
  OPMisc::getMeankT() const
  {
    return 2.0 * _KE.mean() / (Sim->N * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getMeanSqrkT() const
  {
    return 4.0 * _KE.meanSqr()
      / (Sim->N * Sim->N * Sim->dynamics->getParticleDOF()
	 * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getCurrentkT() const
  {
    return 2.0 * _KE.current() / (Sim->N * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getMeanUConfigurational() const
  { 
    return _internalE.mean(); 
  }

  double 
  OPMisc::getMeanSqrUConfigurational() const
  { return _internalE.meanSqr(); }

  void
  OPMisc::initialise()
  {
    _KE.init(Sim->dynamics->getSystemKineticEnergy());
    _internalE.init(Sim->calcInternalEnergy());

    dout << "Particle Count " << Sim->N
	 << "\nSim Unit Length " << Sim->units.unitLength()
	 << "\nSim Unit Time " << Sim->units.unitTime()
	 << "\nDensity " << Sim->getNumberDensity()
      * Sim->units.unitVolume()
	 << "\nPacking Fraction " << Sim->getPackingFraction()
	 << "\nTemperature " << getCurrentkT() / Sim->units.unitEnergy() << std::endl;

    dout << "No. of Species " << Sim->species.size()
	 << "\nSimulation box length <x y z> < ";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << Sim->primaryCellSize[iDim] / Sim->units.unitLength() << " ";
    dout << ">" << std::endl;

    collisionalP.zero();

    Matrix kineticPinit;
    kineticPinit.zero();

    Vector sysMomentum(0,0,0);
    Vector thermalConductivityFS(0,0,0);
    BOOST_FOREACH(const Particle& part, Sim->particles)
      {
	double mass = Sim->species[part]->getMass(part.getID());
	kineticPinit += mass * Dyadic(part.getVelocity(),part.getVelocity());
	sysMomentum += mass * part.getVelocity();
	thermalConductivityFS += part.getVelocity () * Sim->dynamics->getParticleKineticEnergy(part);
      }

    _kineticP.init(kineticPinit);
    _sysMomentum.init(sysMomentum);

    //Set up the correlators
    double correlator_dt = Sim->lastRunMFT / 8;
    if (correlator_dt == 0.0)
      correlator_dt = 1.0 / sqrt(getCurrentkT());
    
    _thermalConductivity.resize(correlator_dt, 10);
    _thermalConductivity.setFreeStreamValue(thermalConductivityFS);

    dout << "Total momentum <x,y,z> < ";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << sysMomentum[iDim] / Sim->units.unitMomentum() << " ";
    dout << ">" << std::endl;

    std::time(&tstartTime);

    clock_gettime(CLOCK_MONOTONIC, &acc_tstartTime);

    std::string sTime(std::ctime(&tstartTime));
    sTime[sTime.size()-1] = ' ';

    dout << "Started on " << sTime << std::endl;
  }

  void
  OPMisc::eventUpdate(const IntEvent& eevent, const PairEventData& PDat)
  {
    stream(eevent.getdt());
    eventUpdate(PDat);
  }

  void
  OPMisc::eventUpdate(const GlobalEvent& eevent, const NEventData& NDat)
  {
    stream(eevent.getdt());
    eventUpdate(NDat);
  }

  void
  OPMisc::eventUpdate(const LocalEvent& eevent, const NEventData& NDat)
  {
    stream(eevent.getdt());
    eventUpdate(NDat);
  }

  void
  OPMisc::eventUpdate(const System&, const NEventData& NDat,
		      const double& dt)
  {
    stream(dt);
    eventUpdate(NDat);
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
  }

  void OPMisc::eventUpdate(const NEventData& NDat)
  {
    BOOST_FOREACH(const ParticleEventData& PDat, NDat.L1partChanges)
      {
        _singleEvents += (PDat.getType() != VIRTUAL);
	_virtualEvents += (PDat.getType() == VIRTUAL);
	const Particle& part = Sim->particles[PDat.getParticleID()];
	const double p1E = Sim->dynamics->getParticleKineticEnergy(Sim->particles[PDat.getParticleID()]);
	
	_KE += PDat.getDeltaKE();
	_internalE += PDat.getDeltaU();
	
	double mass = Sim->species[PDat.getSpeciesID()]->getMass(part.getID());
	_kineticP
	  += mass
	  * (Dyadic(part.getVelocity(), part.getVelocity())
	     - Dyadic(PDat.getOldVel(), PDat.getOldVel()));

	_sysMomentum += mass * (part.getVelocity() - PDat.getOldVel());

	_thermalConductivity.setFreeStreamValue
	  (_thermalConductivity.getFreeStreamValue() 
	   + part.getVelocity() * p1E
	   - PDat.getOldVel() * (p1E - PDat.getDeltaKE()));
      }

    BOOST_FOREACH(const PairEventData& PDat, NDat.L2partChanges)
      eventUpdate(PDat);
  }

  void OPMisc::eventUpdate(const PairEventData& PDat)
  {
    _dualEvents += (PDat.getType() != VIRTUAL);
    _virtualEvents += (PDat.getType() == VIRTUAL);

    _KE += PDat.particle1_.getDeltaKE() + PDat.particle2_.getDeltaKE();
    _internalE += PDat.particle1_.getDeltaU() + PDat.particle2_.getDeltaU();

    const Particle& part1 = Sim->particles[PDat.particle1_.getParticleID()];
    const Particle& part2 = Sim->particles[PDat.particle2_.getParticleID()];
    const Species& sp1 = *Sim->species[PDat.particle1_.getSpeciesID()];
    const Species& sp2 = *Sim->species[PDat.particle2_.getSpeciesID()];
    const double p1E = Sim->dynamics->getParticleKineticEnergy(part1);
    const double p2E = Sim->dynamics->getParticleKineticEnergy(part2);
    
    Vector delP = sp1.getMass(part1.getID()) * (part1.getVelocity() - PDat.particle1_.getOldVel());

    collisionalP += Dyadic(delP, PDat.rij);

    _kineticP
      += sp1.getMass(part1.getID())
      * (Dyadic(part1.getVelocity(), part1.getVelocity())
	 - Dyadic(PDat.particle1_.getOldVel(), PDat.particle1_.getOldVel()))
      + sp2.getMass(part2.getID())
      * (Dyadic(part2.getVelocity(), part2.getVelocity())
	 - Dyadic(PDat.particle2_.getOldVel(), PDat.particle2_.getOldVel()));

    _thermalConductivity.addImpulse(PDat.rij * PDat.particle1_.getDeltaKE());

    _thermalConductivity.setFreeStreamValue
      (_thermalConductivity.getFreeStreamValue() 
       + part1.getVelocity() * p1E + part2.getVelocity() * p2E
       - PDat.particle1_.getOldVel() * (p1E - PDat.particle1_.getDeltaKE())
       - PDat.particle2_.getOldVel() * (p2E - PDat.particle2_.getDeltaKE()));
  }

  double
  OPMisc::getMFT() const
  {
    return Sim->systemTime * static_cast<double>(Sim->N)
      /(Sim->units.unitTime()
	* ((2.0 * static_cast<double>(_dualEvents))
	   + static_cast<double>(_singleEvents)));
  }

  double 
  OPMisc::getEventsPerSecond() const
  {
    timespec acc_tendTime;
    clock_gettime(CLOCK_MONOTONIC, &acc_tendTime);
    
    double duration = double(acc_tendTime.tv_sec) - double(acc_tstartTime.tv_sec)
      + 1e-9 * (double(acc_tendTime.tv_nsec) - double(acc_tstartTime.tv_nsec));

    return Sim->eventCount / duration;
  }

  double 
  OPMisc::getSimTimePerSecond() const
  {
    timespec acc_tendTime;
    clock_gettime(CLOCK_MONOTONIC, &acc_tendTime);
    
    double duration = double(acc_tendTime.tv_sec) - double(acc_tstartTime.tv_sec)
      + 1e-9 * (double(acc_tendTime.tv_nsec) - double(acc_tstartTime.tv_nsec));

    return Sim->systemTime / (duration * Sim->units.unitTime());
  }


  void
  OPMisc::output(magnet::xml::XmlStream &XML)
  {
    std::time_t tendTime;
    time(&tendTime);

    std::string sTime(std::ctime(&tstartTime));
    //A hack to remove the newline character at the end
    sTime[sTime.size()-1] = ' ';

    std::string eTime(std::ctime(&tendTime));
    //A hack to remove the newline character at the end
    eTime[eTime.size()-1] = ' ';

    dout << "Ended on " << eTime
	 << "\nTotal Collisions Executed " << Sim->eventCount
	 << "\nAvg Events/s " << getEventsPerSecond()
	 << "\nSim time per second " << getSimTimePerSecond()
	 << std::endl;

    XML << magnet::xml::tag("Misc")
	<< magnet::xml::tag("Density")
	<< magnet::xml::attr("val")
	<< Sim->getNumberDensity() * Sim->units.unitVolume()
	<< magnet::xml::endtag("Density")

	<< magnet::xml::tag("PackingFraction")
	<< magnet::xml::attr("val") << Sim->getPackingFraction()
	<< magnet::xml::endtag("PackingFraction")

	<< magnet::xml::tag("SpeciesCount")
	<< magnet::xml::attr("val") << Sim->species.size()
	<< magnet::xml::endtag("SpeciesCount")

	<< magnet::xml::tag("ParticleCount")
	<< magnet::xml::attr("val") << Sim->N
	<< magnet::xml::endtag("ParticleCount")

	<< magnet::xml::tag("SystemMomentum")
	<< magnet::xml::tag("Current")
	<< magnet::xml::attr("x") << _sysMomentum.current()[0] / Sim->units.unitMomentum()
	<< magnet::xml::attr("y") << _sysMomentum.current()[1] / Sim->units.unitMomentum()
	<< magnet::xml::attr("z") << _sysMomentum.current()[2] / Sim->units.unitMomentum()
	<< magnet::xml::endtag("Current")
	<< magnet::xml::tag("Average")
	<< magnet::xml::attr("x") << _sysMomentum.mean()[0] / Sim->units.unitMomentum()
	<< magnet::xml::attr("y") << _sysMomentum.mean()[1] / Sim->units.unitMomentum()
	<< magnet::xml::attr("z") << _sysMomentum.mean()[2] / Sim->units.unitMomentum()
	<< magnet::xml::endtag("Average")
	<< magnet::xml::endtag("SystemMomentum")

	<< magnet::xml::tag("Temperature")
	<< magnet::xml::attr("Mean") << getMeankT() / Sim->units.unitEnergy()
	<< magnet::xml::attr("MeanSqr") << getMeanSqrkT() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Current") << getCurrentkT() / Sim->units.unitEnergy()
	<< magnet::xml::attr("Min") << 2.0 * _KE.min() / (Sim->N * Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Max") << 2.0 * _KE.max() / (Sim->N * Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	<< magnet::xml::endtag("Temperature")

	<< magnet::xml::tag("UConfigurational")
	<< magnet::xml::attr("Mean") << getMeanUConfigurational() / Sim->units.unitEnergy()
	<< magnet::xml::attr("MeanSqr") << getMeanSqrUConfigurational() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Current") << _internalE.current() / Sim->units.unitEnergy()
	<< magnet::xml::attr("Min") << _internalE.min() / Sim->units.unitEnergy()
	<< magnet::xml::attr("Max") << _internalE.max() / Sim->units.unitEnergy()
	<< magnet::xml::endtag("UConfigurational")

	<< magnet::xml::tag("ResidualHeatCapacity")
	<< magnet::xml::attr("Value") 
	<< (getMeanSqrUConfigurational() - getMeanUConfigurational() * getMeanUConfigurational())
      / (getMeankT() * getMeankT())
	<< magnet::xml::endtag("ResidualHeatCapacity")
	<< magnet::xml::tag("Pressure")
	<< magnet::xml::attr("Avg") << (_kineticP.mean().tr() + collisionalP.tr() / Sim->systemTime)
	    / (3.0 * Sim->getSimVolume() / Sim->units.unitPressure())
	<< magnet::xml::tag("Tensor") << magnet::xml::chardata()
      ;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  XML << (_kineticP.mean()(iDim, jDim) + collisionalP(iDim, jDim) / Sim->systemTime)
	    / (Sim->getSimVolume() / Sim->units.unitPressure())
	      << " ";
	XML << "\n";
      }
    
    XML << magnet::xml::endtag("Tensor")
	<< magnet::xml::endtag("Pressure")
	<< magnet::xml::tag("Duration")
	<< magnet::xml::attr("Events") << Sim->eventCount
	<< magnet::xml::attr("OneParticleEvents") << _singleEvents
	<< magnet::xml::attr("TwoParticleEvents") << _dualEvents
	<< magnet::xml::attr("VirtualEvents") << _virtualEvents
	<< magnet::xml::attr("Time") << Sim->systemTime / Sim->units.unitTime()
	<< magnet::xml::endtag("Duration")

	<< magnet::xml::tag("Timing")
	<< magnet::xml::attr("Start") << sTime
	<< magnet::xml::attr("End") << eTime
	<< magnet::xml::attr("EventsPerSec") << getEventsPerSecond()
	<< magnet::xml::attr("SimTimePerSec") << getSimTimePerSecond()
	<< magnet::xml::endtag("Timing")

	<< magnet::xml::tag("PrimaryImageSimulationSize")
	<< Sim->primaryCellSize / Sim->units.unitLength()
	<< magnet::xml::endtag("PrimaryImageSimulationSize");

    Vector sumMV(0, 0, 0);
    //Determine the system momentum
    BOOST_FOREACH( const Particle & Part, Sim->particles)
      sumMV += Part.getVelocity() * Sim->species[Part]->getMass(Part.getID());

    XML << magnet::xml::tag("Total_momentum")
	<< sumMV / Sim->units.unitMomentum()
	<< magnet::xml::endtag("Total_momentum")
	<< magnet::xml::tag("totMeanFreeTime")
	<< magnet::xml::attr("val")
	<< getMFT()
	<< magnet::xml::endtag("totMeanFreeTime")
	<< magnet::xml::tag("NegativeTimeEvents")
	<< magnet::xml::attr("Count") << _reverseEvents
	<< magnet::xml::endtag("NegativeTimeEvents")
	<< magnet::xml::tag("Memusage")
	<< magnet::xml::attr("MaxKiloBytes") << magnet::process_mem_usage()
	<< magnet::xml::endtag("Memusage")
	<< magnet::xml::tag("ThermalConductivityCorrelator")
	<< magnet::xml::chardata();

    std::vector<magnet::math::LogarithmicTimeCorrelator<Vector>::Data>
      thermaldata = _thermalConductivity.getAveragedCorrelator();
    
    double inv_units = Sim->units.unitk()
      / ( Sim->units.unitTime() * Sim->units.unitThermalCond() * 2.0 
	  * std::pow(getMeankT(), 2) * Sim->getSimVolume());

    XML << "0 0 0 0 0\n";
    for (size_t i(0); i < thermaldata.size(); ++i)
      XML << thermaldata[i].time / Sim->units.unitTime() << " "
	  << thermaldata[i].sample_count << " "
	  << thermaldata[i].value[0] * inv_units << " "
	  << thermaldata[i].value[1] * inv_units << " "
	  << thermaldata[i].value[2] * inv_units << "\n";

    XML << magnet::xml::endtag("ThermalConductivityCorrelator")
	<< magnet::xml::endtag("Misc");
  }

  void
  OPMisc::periodicOutput()
  {
    time_t rawtime;
    time(&rawtime);

    tm timeInfo;
    localtime_r (&rawtime, &timeInfo);

    char dateString[12] = "";
    strftime(dateString, 12, "%a %H:%M", &timeInfo);

    //Output the date
    I_Pcout() << dateString;

    //Calculate the ETA of the simulation, and take care with overflows and the like
    double _earliest_end_time = HUGE_VAL;
    BOOST_FOREACH(const std::tr1::shared_ptr<System>& sysPtr, Sim->systems)
      if (std::tr1::dynamic_pointer_cast<SystHalt>(sysPtr))
	_earliest_end_time = std::min(_earliest_end_time, sysPtr->getdt());

    double time_seconds_remaining = _earliest_end_time 
      / (getSimTimePerSecond() * Sim->units.unitTime());

    size_t seconds_remaining = time_seconds_remaining;
    
    if (time_seconds_remaining > std::numeric_limits<size_t>::max())
      seconds_remaining = std::numeric_limits<size_t>::max();

    if (Sim->endEventCount != std::numeric_limits<unsigned long long>::max())
      {
	double event_seconds_remaining = (Sim->endEventCount - Sim->eventCount) / getEventsPerSecond() + 0.5;
	
	if (event_seconds_remaining < std::numeric_limits<size_t>::max())
	  seconds_remaining = std::min(seconds_remaining, size_t(event_seconds_remaining));

      }
  
    if (seconds_remaining != std::numeric_limits<size_t>::max())
      {
	size_t ETA_hours = seconds_remaining / 3600;
	size_t ETA_mins = (seconds_remaining / 60) % 60;
	size_t ETA_secs = seconds_remaining % 60;

	I_Pcout() << ", ETA ";
	if (ETA_hours)
	  I_Pcout() << ETA_hours << "hr ";
	  
	if (ETA_mins)
	  I_Pcout() << ETA_mins << "min ";

	I_Pcout() << ETA_secs << "s";
      }

    I_Pcout() << ", Events " << (Sim->eventCount+1)/1000 << "k, t "
	      << Sim->systemTime/Sim->units.unitTime() 
	      << ", <Mean Free Time> " <<  getMFT()
	      << ", ";
  }
}

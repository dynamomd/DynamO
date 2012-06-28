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
#include <boost/foreach.hpp>
#include <sys/time.h>
#include <ctime>

namespace dynamo {
  OPMisc::OPMisc(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"Misc",0),
    oldSysTime(0),
    dualEvents(0),
    singleEvents(0),
    oldcoll(0),
    _reverseEvents(0),
    InitialKE(0.0),
    KEacc(0.0),
    KEsqAcc(0.0),
    KECurrent(0.0),
    intECurrent(0.0),
    intEsqAcc(0.0),
    intEAcc(0.0)
  {}

  void
  OPMisc::changeSystem(OutputPlugin* misc2)
  {
    OPMisc& op = static_cast<OPMisc&>(*misc2);
    std::swap(Sim, op.Sim);
    std::swap(KECurrent, op.KECurrent);
    std::swap(intECurrent, op.intECurrent);
    std::swap(curr_kineticP, op.curr_kineticP);
  }

  void
  OPMisc::temperatureRescale(const double& scale)
  { KECurrent *= scale; }

  double 
  OPMisc::getMeankT() const
  {
    if (KEacc == 0) return getCurrentkT();

    return 2.0 * KEacc / (Sim->dSysTime * Sim->N * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getMeanSqrkT() const
  {
    if (KEsqAcc == 0) 
      return 4.0 * KECurrent * KECurrent
	/ (Sim->N * Sim->N
	   * Sim->dynamics->getParticleDOF()
	   * Sim->dynamics->getParticleDOF());
    
    return 4.0 * KEsqAcc 
	/ (Sim->dSysTime
	   * Sim->N * Sim->N
	   * Sim->dynamics->getParticleDOF()
	   * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getCurrentkT() const
  {
    return 2.0 * KECurrent / (Sim->N * Sim->dynamics->getParticleDOF());
  }

  double 
  OPMisc::getMeanUConfigurational() const
  { 
    return intEAcc / Sim->dSysTime; 
  }

  double 
  OPMisc::getMeanSqrUConfigurational() const
  { return intEsqAcc / Sim->dSysTime; }

  void
  OPMisc::initialise()
  {
    KEMax = KEMin = InitialKE = KECurrent = Sim->dynamics->getSystemKineticEnergy();
    intEMax = intEMin = intECurrent = Sim->calcInternalEnergy();

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
      dout  << Sim->primaryCellSize[iDim]/Sim->units.unitLength() << " ";
    dout << ">" << std::endl;

    Vector sumMV (0,0,0);

    //Determine the discrepancy VECTOR
    BOOST_FOREACH( const Particle & Part, Sim->particleList)
      {
	Vector  pos(Part.getPosition()), vel(Part.getVelocity());
	Sim->BCs->applyBC(pos, vel);
	sumMV += vel * Sim->species[Part]->getMass(Part.getID());
      }

    dout << "Total momentum <x,y,z> <";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << sumMV[iDim] / Sim->units.unitMomentum() << " ";
    dout << ">" << std::endl;

    cumulative_kineticP.zero();
    collisionalP.zero();
    curr_kineticP.zero();
    
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	curr_kineticP 
	  += Sim->species[part]->getMass(part.getID())
	  * Dyadic(part.getVelocity(),part.getVelocity());
      }

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
    KEacc += KECurrent * dt;
    KEsqAcc += KECurrent * KECurrent * dt;
    intEAcc += intECurrent * dt;
    intEsqAcc += intECurrent * intECurrent * dt;
    cumulative_kineticP += curr_kineticP * dt;
  }

  void OPMisc::eventUpdate(const NEventData& NDat)
  {
    BOOST_FOREACH(const ParticleEventData& PDat, NDat.L1partChanges)
      {
	++singleEvents;
	const Particle& part = Sim->particleList[PDat.getParticleID()];
	
	KECurrent += PDat.getDeltaKE();
	intECurrent += PDat.getDeltaU();
	
	curr_kineticP
	  += Sim->species[PDat.getSpeciesID()]->getMass(part.getID())
	  * (Dyadic(part.getVelocity(), part.getVelocity())
	     - Dyadic(PDat.getOldVel(), PDat.getOldVel()));
      }

    BOOST_FOREACH(const PairEventData& PDat, NDat.L2partChanges)
      eventUpdate(PDat);

    KEMin = std::min(KECurrent, KEMin);
    KEMax = std::max(KECurrent, KEMax);
    intEMin = std::min(intECurrent, intEMin);
    intEMax = std::max(intECurrent, intEMax);
  }

  void OPMisc::eventUpdate(const PairEventData& PDat)
  {
    if ((PDat.getType() != NBHOOD_IN)
	&& (PDat.getType() != NBHOOD_OUT))
      ++dualEvents;

    KECurrent += PDat.particle1_.getDeltaKE() + PDat.particle2_.getDeltaKE();
    intECurrent += PDat.particle1_.getDeltaU() + PDat.particle2_.getDeltaU();

    const Particle& part1 = Sim->particleList[PDat.particle1_.getParticleID()];
    const Particle& part2 = Sim->particleList[PDat.particle2_.getParticleID()];
    const Species& sp1 = *Sim->species[PDat.particle1_.getSpeciesID()];
    const Species& sp2 = *Sim->species[PDat.particle2_.getSpeciesID()];
    
    Vector delP = sp1.getMass(part1.getID()) * (part1.getVelocity() - PDat.particle1_.getOldVel());

    collisionalP += Dyadic(delP, PDat.rij);

    curr_kineticP 
      += sp1.getMass(part1.getID())
      * (Dyadic(part1.getVelocity(), part1.getVelocity())
	 - Dyadic(PDat.particle1_.getOldVel(), PDat.particle1_.getOldVel()))
      + sp2.getMass(part2.getID())
      * (Dyadic(part2.getVelocity(), part2.getVelocity())
	 - Dyadic(PDat.particle2_.getOldVel(), PDat.particle2_.getOldVel()))
      ;

    KEMin = std::min(KECurrent, KEMin);
    KEMax = std::max(KECurrent, KEMax);
    intEMin = std::min(intECurrent, intEMin);
    intEMax = std::max(intECurrent, intEMax);
  }

  double
  OPMisc::getMFT() const
  {
    return Sim->dSysTime * static_cast<double>(Sim->N)
      /(Sim->units.unitTime()
	* ((2.0 * static_cast<double>(dualEvents))
	   + static_cast<double>(singleEvents)));
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

    return Sim->dSysTime / (duration * Sim->units.unitTime());
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

	<< magnet::xml::tag("Temperature")
	<< magnet::xml::attr("Mean") << getMeankT() / Sim->units.unitEnergy()
	<< magnet::xml::attr("MeanSqr") << getMeanSqrkT() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Current") << getCurrentkT() / Sim->units.unitEnergy()
	<< magnet::xml::attr("Min") << 2.0 * KEMin / (Sim->N * Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Max") << 2.0 * KEMax / (Sim->N * Sim->dynamics->getParticleDOF() * Sim->units.unitEnergy())
	<< magnet::xml::endtag("Temperature")

	<< magnet::xml::tag("UConfigurational")
	<< magnet::xml::attr("Mean") << getMeanUConfigurational() / Sim->units.unitEnergy()
	<< magnet::xml::attr("MeanSqr") << getMeanSqrUConfigurational() / (Sim->units.unitEnergy() * Sim->units.unitEnergy())
	<< magnet::xml::attr("Current") << intECurrent / Sim->units.unitEnergy()
	<< magnet::xml::attr("Min") << intEMin / Sim->units.unitEnergy()
	<< magnet::xml::attr("Max") << intEMax / Sim->units.unitEnergy()
	<< magnet::xml::endtag("UConfigurational")

	<< magnet::xml::tag("ResidualHeatCapacity")
	<< magnet::xml::attr("Value") 
	<< (getMeanSqrUConfigurational() - getMeanUConfigurational() * getMeanUConfigurational())
      / (getMeankT() * getMeankT())
	<< magnet::xml::endtag("ResidualHeatCapacity")
	<< magnet::xml::tag("Pressure")
	<< magnet::xml::attr("Avg") << (cumulative_kineticP.tr() + collisionalP.tr())
	    / (3.0 * Sim->dSysTime * Sim->getSimVolume() / Sim->units.unitPressure())
	<< magnet::xml::tag("Tensor") << magnet::xml::chardata()
      ;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  XML << (cumulative_kineticP(iDim, jDim) + collisionalP(iDim, jDim))
	    / (Sim->dSysTime * Sim->getSimVolume() / Sim->units.unitPressure())
	      << " ";
	XML << "\n";
      }
    
    XML << magnet::xml::endtag("Tensor")
	<< magnet::xml::endtag("Pressure")
	<< magnet::xml::tag("Duration")
	<< magnet::xml::attr("Events") << Sim->eventCount
	<< magnet::xml::attr("OneParticleEvents") << singleEvents
	<< magnet::xml::attr("TwoParticleEvents") << dualEvents
	<< magnet::xml::attr("Time") << Sim->dSysTime / Sim->units.unitTime()
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
    BOOST_FOREACH( const Particle & Part, Sim->particleList)
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

    //Output when the simulation will end
    if (Sim->endEventCount != std::numeric_limits<unsigned long long>::max())
      {
	size_t seconds_remaining = (Sim->endEventCount - Sim->eventCount) / getEventsPerSecond() + 0.5;
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
	      << Sim->dSysTime/Sim->units.unitTime() 
	      << ", <Mean Free Time> " <<  getMFT()
	      << ", ";

    oldSysTime = Sim->dSysTime;
    oldcoll = Sim->eventCount;
  }
}

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

#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <magnet/memUsage.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>
#include <sys/time.h>
#include <ctime>

namespace dynamo {
  OPMisc::OPMisc(const dynamo::SimData* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"Misc",0),
    oldSysTime(0),
    dualEvents(0),
    singleEvents(0),
    oldcoll(0)
  {}

  void
  OPMisc::changeSystem(OutputPlugin* misc2)
  {
    std::swap(Sim, static_cast<OPMisc*>(misc2)->Sim);
  }

  void
  OPMisc::initialise()
  {
    double kt = Sim->dynamics.getLiouvillean().getkT();


    dout << "Particle Count " << Sim->N
	 << "\nSim Unit Length " << Sim->dynamics.units().unitLength()
	 << "\nSim Unit Time " << Sim->dynamics.units().unitTime()
	 << "\nDensity " << Sim->dynamics.getNumberDensity()
      * Sim->dynamics.units().unitVolume()
	 << "\nPacking Fraction " << Sim->dynamics.getPackingFraction()
	 << "\nSim Temperature " << kt
	 << "\nReduced Temperature " << kt / Sim->dynamics.units().unitEnergy() << std::endl;

    dout << "No. of Species " << Sim->dynamics.getSpecies().size()
	 << "\nSimulation box length <x y z> < ";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << Sim->primaryCellSize[iDim]/Sim->dynamics.units().unitLength() << " ";
    dout << ">" << std::endl;

    Vector sumMV (0,0,0);

    //Determine the discrepancy VECTOR
    BOOST_FOREACH( const Particle & Part, Sim->particleList)
      {
	Vector  pos(Part.getPosition()), vel(Part.getVelocity());
	Sim->dynamics.BCs().applyBC(pos, vel);
	sumMV += vel * Sim->dynamics.getSpecies(Part).getMass(Part.getID());
      }

    dout << "Total momentum <x,y,z> <";
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      dout  << sumMV[iDim] / Sim->dynamics.units().unitMomentum() << " ";
    dout << ">" << std::endl;

    std::time(&tstartTime);

    clock_gettime(CLOCK_MONOTONIC, &acc_tstartTime);

    std::string sTime(std::ctime(&tstartTime));
    sTime[sTime.size()-1] = ' ';

    dout << "Started on " << sTime << std::endl;
  }

  void
  OPMisc::eventUpdate(const IntEvent&, const PairEventData&)
  {
    ++dualEvents;
  }

  void
  OPMisc::eventUpdate(const GlobalEvent&, const NEventData& NDat)
  {
    dualEvents += NDat.L2partChanges.size();
    singleEvents += NDat.L1partChanges.size();
  }

  void
  OPMisc::eventUpdate(const LocalEvent&, const NEventData& NDat)
  {
    dualEvents += NDat.L2partChanges.size();
    singleEvents += NDat.L1partChanges.size();
  }

  void
  OPMisc::eventUpdate(const System&, const NEventData& NDat,
		      const double&)
  {
    dualEvents += NDat.L2partChanges.size();
    singleEvents += NDat.L1partChanges.size();
  }

  double
  OPMisc::getMFT() const
  {
    return Sim->dSysTime * static_cast<double>(Sim->N)
      /(Sim->dynamics.units().unitTime()
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

    return Sim->dSysTime / (duration * Sim->dynamics.units().unitTime());
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
	<< magnet::xml::tag("Memusage")
	<< magnet::xml::attr("MaxKiloBytes") << magnet::process_mem_usage()
	<< magnet::xml::endtag("Memusage")
	<< magnet::xml::tag("Density")
	<< magnet::xml::attr("val")
	<< Sim->dynamics.getNumberDensity() * Sim->dynamics.units().unitVolume()
	<< magnet::xml::endtag("Density")

	<< magnet::xml::tag("PackingFraction")
	<< magnet::xml::attr("val") << Sim->dynamics.getPackingFraction()
	<< magnet::xml::endtag("PackingFraction")

	<< magnet::xml::tag("SpeciesCount")
	<< magnet::xml::attr("val") << Sim->dynamics.getSpecies().size()
	<< magnet::xml::endtag("SpeciesCount")

	<< magnet::xml::tag("ParticleCount")
	<< magnet::xml::attr("val") << Sim->N
	<< magnet::xml::endtag("ParticleCount")

	<< magnet::xml::tag("Duration")
	<< magnet::xml::attr("Events") << Sim->eventCount
	<< magnet::xml::attr("OneParticleEvents") << singleEvents
	<< magnet::xml::attr("TwoParticleEvents") << dualEvents
	<< magnet::xml::attr("Time") << Sim->dSysTime / Sim->dynamics.units().unitTime()
	<< magnet::xml::endtag("Duration")

	<< magnet::xml::tag("Timing")
	<< magnet::xml::attr("Start") << sTime
	<< magnet::xml::attr("End") << eTime
	<< magnet::xml::attr("EventsPerSec") << getEventsPerSecond()
	<< magnet::xml::attr("SimTimePerSec") << getSimTimePerSecond()
	<< magnet::xml::endtag("Timing")

	<< magnet::xml::tag("SystemBoxLength")
	<< Sim->primaryCellSize / Sim->dynamics.units().unitLength()
	<< magnet::xml::endtag("SystemBoxLength");

    Vector sumMV(0, 0, 0);
    //Determine the system momentum
    BOOST_FOREACH( const Particle & Part, Sim->particleList)
      sumMV += Part.getVelocity() * Sim->dynamics.getSpecies(Part).getMass(Part.getID());

    XML << magnet::xml::tag("Total_momentum")
	<< sumMV / Sim->dynamics.units().unitMomentum()
	<< magnet::xml::endtag("Total_momentum")
	<< magnet::xml::tag("totMeanFreeTime")
	<< magnet::xml::attr("val")
	<< getMFT()
	<< magnet::xml::endtag("totMeanFreeTime")
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
	      << Sim->dSysTime/Sim->dynamics.units().unitTime() << ", <t_2> "
	      <<   Sim->dSysTime * static_cast<double>(Sim->N)
      /(Sim->dynamics.units().unitTime() * 2.0 * static_cast<double>(dualEvents))
	      << ", <t_tot> "
	      <<   Sim->dSysTime * static_cast<double>(Sim->N)
      / (Sim->dynamics.units().unitTime() * (2.0 * static_cast<double>(dualEvents)
					     + static_cast<double>(singleEvents)))
	      << ", ";

    oldSysTime = Sim->dSysTime;
    oldcoll = Sim->eventCount;
  }
}

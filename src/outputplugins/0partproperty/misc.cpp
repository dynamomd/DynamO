/*  DYNAMO:- Event driven molecular dynamics simulator
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

#include "misc.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../datatypes/vector.xml.hpp"
#include <ctime>
#include <sys/time.h>
#include <magnet/memUsage.hpp>
#include "../../dynamics/species/inertia.hpp"

OPMisc::OPMisc(const DYNAMO::SimData* tmp, const XMLNode&):
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

  Vector  VecEnergy(Sim->dynamics.getLiouvillean().getVectorSystemKineticEnergy());

  VecEnergy *= 2.0 / (Sim->N * Sim->dynamics.units().unitEnergy());

  I_cout() << "Particle Count " << Sim->N
	   << "\nSim Unit Length " << Sim->dynamics.units().unitLength()
	   << "\nSim Unit Time " << Sim->dynamics.units().unitTime()
	   << "\nDensity " << Sim->dynamics.getNumberDensity()
    * Sim->dynamics.units().unitVolume()
	   << "\nPacking Fraction " << Sim->dynamics.getPackingFraction()
	   << "\nSim Temperature " << kt
	   << "\nReduced Temperature " << kt / Sim->dynamics.units().unitEnergy();

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    I_cout() << "Kinetic Temperature dimension" << iDim << " "
	     <<  VecEnergy[iDim];

  I_cout() << "No. of Species " << Sim->dynamics.getSpecies().size()
      << "\nSimulation box length <x,y,z> ";
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    std::cout  << Sim->aspectRatio[iDim]/Sim->dynamics.units().unitLength() << " ";

  Vector  sumMV (0,0,0);

  //Determine the discrepancy VECTOR
  BOOST_FOREACH( const Particle & Part, Sim->particleList)
    {
      Vector  pos(Part.getPosition()), vel(Part.getVelocity());
      Sim->dynamics.BCs().applyBC(pos, vel);

      sumMV += vel * Sim->dynamics.getSpecies(Part).getMass();
    }

  I_cout() << "Total momentum <x,y,z> <";

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    std::cout  << sumMV[iDim] / Sim->dynamics.units().unitMomentum() << " ";

  std::cout << ">";

  std::time(&tstartTime);

  clock_gettime(CLOCK_MONOTONIC, &acc_tstartTime);

  std::string sTime(std::ctime(&tstartTime));
  sTime[sTime.size()-1] = ' ';

  I_cout() << "Started on " << sTime;
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


void
OPMisc::output(xml::XmlStream &XML)
{
  std::time_t tendTime;
  time(&tendTime);

  std::string sTime(std::ctime(&tstartTime));
  //A hack to remove the newline character at the end
  sTime[sTime.size()-1] = ' ';

  std::string eTime(std::ctime(&tendTime));
  //A hack to remove the newline character at the end
  eTime[eTime.size()-1] = ' ';

  timespec acc_tendTime;
  clock_gettime(CLOCK_MONOTONIC, &acc_tendTime);

  double collpersec = static_cast<double>(Sim->eventCount)
    / (double(acc_tendTime.tv_sec) + 1e-9 * double(acc_tendTime.tv_nsec)
       - double(acc_tstartTime.tv_sec) - 1e-9 * double(acc_tstartTime.tv_nsec));

  I_cout() << "Ended on " << eTime
	   << "\nTotal Collisions Executed " << Sim->eventCount
	   << "\nAvg Coll/s " << collpersec
	   << "\nSim time per second "
	   << Sim->dSysTime / (Sim->dynamics.units().unitTime()
			       * static_cast<double>(tendTime - tstartTime));

  XML << xml::tag("Misc")
      << xml::tag("Memusage")
      << xml::attr("MaxKiloBytes") << magnet::process_mem_usage()
      << xml::endtag("Memusage")
      << xml::tag("Density")
      << xml::attr("val") << Sim->dynamics.getNumberDensity() * Sim->dynamics.units().unitVolume()
      << xml::endtag("Density")

      << xml::tag("PackingFraction")
      << xml::attr("val") << Sim->dynamics.getPackingFraction()
      << xml::endtag("PackingFraction")

      << xml::tag("SpeciesCount")
      << xml::attr("val") << Sim->dynamics.getSpecies().size()
      << xml::endtag("SpeciesCount")

      << xml::tag("ParticleCount")
      << xml::attr("val") << Sim->N
      << xml::endtag("ParticleCount")

      << xml::tag("SimLength")
      << xml::attr("Collisions") << Sim->eventCount
      << xml::attr("OneParticleEvents") << singleEvents
      << xml::attr("TwoParticleEvents") << dualEvents
      << xml::attr("Time") << Sim->dSysTime / Sim->dynamics.units().unitTime()
      << xml::endtag("SimLength")

      << xml::tag("Timing")

      << xml::tag("Start")
      << xml::attr("val") << sTime
      << xml::endtag("Start")

      << xml::tag("End")
      << xml::attr("val") << eTime
      << xml::endtag("End")

      << xml::tag("Duration")
      << xml::attr("val")
      << tendTime - tstartTime
      << xml::endtag("Duration")

      << xml::tag("CollPerSec")
      << xml::attr("val") << collpersec
      << xml::attr("CondorWarning") << std::string("true")
      << xml::endtag("CollPerSec")

      << xml::endtag("Timing")
      << xml::tag("SystemBoxLength")
      << xml::attr("val")
      << 1.0/Sim->dynamics.units().unitLength();

  char name[2] = "x";
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      name[0] = 'x' + iDim;
      XML << xml::tag(name) << xml::attr("val")
	  << Sim->aspectRatio[iDim]/Sim->dynamics.units().unitLength()
	  << xml::endtag(name);
    }

  XML << xml::endtag("SystemBoxLength");

  // Output scalar moment of inertia for any species which may have it
  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
  {
    if (dynamic_cast<const SpInertia*>(spec.get_ptr()) != NULL)
    {
      XML << xml::tag("ScalarInertia")
	  << xml::attr("Species") << spec->getName()
	  << xml::attr("Mass") << spec->getMass()
	  << xml::attr("inertiaConst") << (spec->getScalarMomentOfInertia() / (Sim->dynamics.units().unitArea() * spec->getMass()))
	  << xml::endtag("ScalarInertia");
    }
  }

  Vector  sumMV (0,0,0);

  //Determine the discrepancy VECTOR
  BOOST_FOREACH( const Particle & Part, Sim->particleList)
    sumMV += Part.getVelocity() * Sim->dynamics.getSpecies(Part).getMass();

  XML << xml::tag("Total_momentum")
      << sumMV / Sim->dynamics.units().unitMomentum()
      << xml::endtag("Total_momentum")
      << xml::tag("totMeanFreeTime")
      << xml::attr("val")
      << getMFT()
      << xml::endtag("totMeanFreeTime");

  XML << xml::tag("MemoryUsage")
      << xml::attr("ResidentSet") << magnet::process_mem_usage()
      << xml::endtag("MemoryUsage")
      << xml::endtag("Misc");
}

void
OPMisc::periodicOutput()
{
  time_t rawtime;
  time(&rawtime);

  tm timeInfo;
  localtime_r (&rawtime, &timeInfo);

  char dateString[12] = "";
  strftime(dateString, 12, "%a %H:%M |", &timeInfo);

  I_Pcout() << dateString << " NColls " << (Sim->eventCount+1)/1000 << "k, t "
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

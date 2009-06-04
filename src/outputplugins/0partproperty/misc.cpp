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

#include "misc.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../datatypes/vector.xml.hpp"
#include <ctime>

COPMisc::COPMisc(const DYNAMO::SimData* tmp, const XMLNode&):
  COutputPlugin(tmp,"Misc",0),
  oldSysTime(0),
  dualEvents(0),
  singleEvents(0),
  oldcoll(0)
{}

void 
COPMisc::changeSystem(COutputPlugin* misc2)
{
  std::swap(Sim, static_cast<COPMisc*>(misc2)->Sim);
}

void
COPMisc::initialise()
{
  Iflt kt = Sim->Dynamics.Liouvillean().getkT();

  Vector  VecEnergy(Sim->Dynamics.Liouvillean().getVectorSystemKineticEnergy());
  
  VecEnergy *= 2.0 / (Sim->lN * Sim->Dynamics.units().unitEnergy());

  I_cout() << "Particle Count " << Sim->lN 
	   << "\nSim Unit Length " << Sim->Dynamics.units().unitLength()
	   << "\nSim Unit Time " << Sim->Dynamics.units().unitTime()
	   << "\nDensity " << Sim->Dynamics.getNumberDensity() 
    * Sim->Dynamics.units().unitVolume()
	   << "\nPacking Fraction " << Sim->Dynamics.getPackingFraction()
	   << "\nSim Temperature " << kt
	   << "\nReduced Temperature " << kt / Sim->Dynamics.units().unitEnergy();
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    I_cout() << "Kinetic Temperature dimension" << iDim << " " 
	     <<  VecEnergy[iDim];

  I_cout() << "No. of Species " << Sim->Dynamics.getSpecies().size()
      << "\nSimulation box length <x,y,z> ";
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    std::cout  << Sim->aspectRatio[iDim]/Sim->Dynamics.units().unitLength() << " ";
  
  Vector  sumMV (0,0,0);
  
  //Determine the discrepancy VECTOR
  BOOST_FOREACH( const CParticle & Part, Sim->vParticleList)
    {
      Vector  pos(Part.getPosition()), vel(Part.getVelocity());
      Sim->Dynamics.BCs().setPBC(pos, vel);

      sumMV += vel * Sim->Dynamics.getSpecies(Part).getMass();
    }
  
  I_cout() << "Total momentum <x,y,z> <";

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    std::cout  << sumMV[iDim] / Sim->Dynamics.units().unitMomentum() << " ";

  std::cout << ">";

  std::time(&tstartTime);
  
  std::string sTime(std::ctime(&tstartTime));
  sTime[sTime.size()-1] = ' ';

  I_cout() << "Started on " << sTime;
}

void 
COPMisc::eventUpdate(const CIntEvent&, const C2ParticleData&) 
{
  dualEvents++;
}

void 
COPMisc::eventUpdate(const CGlobEvent&, const CNParticleData& NDat)
{
  dualEvents += NDat.L2partChanges.size();
  singleEvents += NDat.L1partChanges.size();
}

void 
COPMisc::eventUpdate(const CLocalEvent&, const CNParticleData& NDat)
{
  dualEvents += NDat.L2partChanges.size();
  singleEvents += NDat.L1partChanges.size();
}

void 
COPMisc::eventUpdate(const CSystem&, const CNParticleData& NDat, 
		     const Iflt&)
{
  dualEvents += NDat.L2partChanges.size();
  singleEvents += NDat.L1partChanges.size();
}

Iflt
COPMisc::getMFT() const
{
  return Sim->dSysTime * static_cast<Iflt>(Sim->lN)
    /(Sim->Dynamics.units().unitTime() 
      * ((2.0 * static_cast<Iflt>(dualEvents)) 
	 + static_cast<Iflt>(singleEvents)));
}


void
COPMisc::output(xmlw::XmlStream &XML)
{
  std::time_t tendTime;
  time(&tendTime);
  
  std::string sTime(std::ctime(&tstartTime));
  //A hack to remove the newline character at the end
  sTime[sTime.size()-1] = ' ';

  std::string eTime(std::ctime(&tendTime));
  //A hack to remove the newline character at the end
  eTime[eTime.size()-1] = ' ';
  
  Iflt collpersec = static_cast<Iflt>(Sim->lNColl)
    / static_cast<Iflt>(tendTime - tstartTime);

  I_cout() << "Ended on " << eTime
	   << "\nTotal Collisions Executed " << Sim->lNColl
	   << "\nAvg Coll/s " << collpersec;
   
  XML << xmlw::tag("Misc")
      << xmlw::tag("Density")
      << xmlw::attr("val") << Sim->Dynamics.getNumberDensity() * Sim->Dynamics.units().unitVolume()
      << xmlw::endtag("Density")
    
      << xmlw::tag("PackingFraction")
      << xmlw::attr("val") << Sim->Dynamics.getPackingFraction()
      << xmlw::endtag("PackingFraction")

      << xmlw::tag("SpeciesCount")
      << xmlw::attr("val") << Sim->Dynamics.getSpecies().size()
      << xmlw::endtag("SpeciesCount")
    
      << xmlw::tag("ParticleCount")
      << xmlw::attr("val") << Sim->lN
      << xmlw::endtag("ParticleCount")
    
      << xmlw::tag("SimLength")
      << xmlw::attr("Collisions") << Sim->lNColl
      << xmlw::attr("Time") << Sim->dSysTime / Sim->Dynamics.units().unitTime()
      << xmlw::endtag("SimLength")
    
      << xmlw::tag("Timing") 

      << xmlw::tag("Start") 
      << xmlw::attr("val") << sTime
      << xmlw::endtag("Start")
  
      << xmlw::tag("End") 
      << xmlw::attr("val") << eTime
      << xmlw::endtag("End")

      << xmlw::tag("Duration") 
      << xmlw::attr("val") 
      << tendTime - tstartTime
      << xmlw::endtag("Duration")

      << xmlw::tag("CollPerSec")
      << xmlw::attr("val") << collpersec
      << xmlw::endtag("CollPerSec")

      << xmlw::endtag("Timing") 
      << xmlw::tag("SystemBoxLength")
      << xmlw::attr("val") 
      << 1.0/Sim->Dynamics.units().unitLength();

  char name[2] = "x";
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      name[0] = 'x' + iDim;
      XML << xmlw::tag(name) << xmlw::attr("val")
	  << Sim->aspectRatio[iDim]/Sim->Dynamics.units().unitLength() 
	  << xmlw::endtag(name);
    }
  
  XML << xmlw::endtag("SystemBoxLength");

  Vector  sumMV (0,0,0);
  
  //Determine the discrepancy VECTOR
  BOOST_FOREACH( const CParticle & Part, Sim->vParticleList)
    sumMV += Part.getVelocity() * Sim->Dynamics.getSpecies(Part).getMass();
  
  XML << xmlw::tag("Total_momentum") 
      << sumMV / Sim->Dynamics.units().unitMomentum()
      << xmlw::endtag("Total_momentum")    
      << xmlw::tag("totMeanFreeTime")
      << xmlw::attr("val")
      << getMFT()
      << xmlw::endtag("totMeanFreeTime")
      << xmlw::endtag("Misc");
}

void
COPMisc::periodicOutput()
{
  time_t rawtime;
  time(&rawtime);

  tm timeInfo;
  localtime_r (&rawtime, &timeInfo);

  char dateString[12] = "";
  strftime(dateString, 12, "%a %H:%M |", &timeInfo);

  I_Pcout() << dateString << " NColls " << (Sim->lNColl+1)/1000 << "k, t "
	    << Sim->dSysTime/Sim->Dynamics.units().unitTime() << ", <t_free> "
	    << getMFT()
	    << ", "; 

  oldSysTime = Sim->dSysTime;
  oldcoll = Sim->lNColl;
}

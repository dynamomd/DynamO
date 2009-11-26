/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "plateMotion.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../base/is_colormap.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/interactions/squarebond.hpp"
#include "../../dynamics/ranges/2RList.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include "../../dynamics/locals/oscillatingplate.hpp"

OPPlateMotion::OPPlateMotion(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"PlateMotion"),
  partpartEnergyLoss(0),
  oldPlateEnergy(0)
{
  operator<<(XML);
}

OPPlateMotion::OPPlateMotion(const OPPlateMotion& cp):
  OPTicker(cp),
  plateID(cp.plateID),
  plateName(cp.plateName),
  partpartEnergyLoss(0),
  oldPlateEnergy(0)
{
  if (cp.logfile.is_open())
    cp.logfile.close();
}

void
OPPlateMotion::initialise()
{
  try {
    plateID = Sim->dynamics.getLocal(plateName)->getID();
  } catch(...)
    {
      D_throw() << "Could not find the PlateName specified. You said " << plateName;
    }
  
  if (dynamic_cast<const CLOscillatingPlate*>(Sim->dynamics.getLocals()[plateID].get_ptr()) == NULL) 
    D_throw() << "The PlateName'd local is not a CLOscillatingPlate";

  if (logfile.is_open())
    logfile.close();
  
  logfile.open("plateMotion.out", std::ios::out|std::ios::trunc);

  localEnergyLoss.resize(Sim->dynamics.getLocals().size(), std::make_pair(Iflt(0),std::vector<Iflt>()));
  localEnergyFlux.resize(Sim->dynamics.getLocals().size(), std::make_pair(Iflt(0),std::vector<Iflt>()));

  oldPlateEnergy = dynamic_cast<const CLOscillatingPlate*>(Sim->dynamics.getLocals()[plateID].get_ptr())->getPlateEnergy();
  partpartEnergyLoss = 0;

  ticker();
}

void 
OPPlateMotion::eventUpdate(const CLocalEvent& localEvent, const CNParticleData& SDat)
{
  Iflt newPlateEnergy = oldPlateEnergy;

  if (localEvent.getLocalID() == plateID)
    newPlateEnergy = dynamic_cast<const CLOscillatingPlate*>(Sim->dynamics.getLocals()[plateID].get_ptr())->getPlateEnergy();

  Iflt EnergyChange(0);

  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    EnergyChange += pData.getDeltaKE();  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    EnergyChange += pData.particle1_.getDeltaKE() + pData.particle2_.getDeltaKE();
  
  localEnergyFlux[localEvent.getLocalID()].first += EnergyChange;

  localEnergyLoss[localEvent.getLocalID()].first += EnergyChange + newPlateEnergy - oldPlateEnergy;

  oldPlateEnergy = newPlateEnergy;
}

void 
OPPlateMotion::eventUpdate(const CIntEvent&, const C2ParticleData& pData)
{
  partpartEnergyLoss += pData.particle1_.getDeltaKE()
    + pData.particle2_.getDeltaKE();
}

void 
OPPlateMotion::ticker()
{
  BOOST_FOREACH(localEntry& entry, localEnergyLoss)
    {
      entry.second.push_back(entry.first);
      entry.first = 0.0;
    }

  Vector com(0,0,0), momentum(0,0,0);
  Iflt sqmom(0);
  Iflt partEnergy(0.0);

  Iflt mass(0);
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      Vector pos(part.getPosition()), vel(part.getVelocity());
      Iflt pmass(Sim->dynamics.getSpecies(part).getMass());
      Sim->dynamics.BCs().applyBC(pos, vel);
      momentum += vel * pmass;
      sqmom += (vel | vel) * (pmass * pmass);
      com += pos * pmass;
      mass += pmass;
      partEnergy += pmass * vel.nrm2();
    }
  
  com /= (mass * Sim->dynamics.units().unitLength());
  Vector comvel = momentum / (mass * Sim->dynamics.units().unitVelocity());
  
  partEnergy *= 0.5;

  const CLOscillatingPlate& plate(*dynamic_cast<const CLOscillatingPlate*>(Sim->dynamics.getLocals()[plateID].get_ptr()));

  Vector platePos = (plate.getPosition() - plate.getCentre()) / Sim->dynamics.units().unitLength();

  Vector plateSpeed = plate.getVelocity() / Sim->dynamics.units().unitVelocity();

  logfile << Sim->dSysTime / Sim->dynamics.units().unitTime()
	  << " " << platePos[0] << " " << platePos[1] << " " << platePos[2] 
	  << " " << com[0] << " " << com[1] << " " << com[2]
	  << " " << comvel[0] << " " << comvel[1] << " " << comvel[2]
	  << " " << plateSpeed[0] << " " << plateSpeed[1] << " " << plateSpeed[2]
	  << " " << (sqmom - ((momentum | momentum) / Sim->lN)) / (Sim->lN * pow(Sim->dynamics.units().unitMomentum(),2))
	  << " " << plate.getPlateEnergy() / Sim->dynamics.units().unitEnergy()
	  << " " << partEnergy / Sim->dynamics.units().unitEnergy() 
	  << " " << (plate.getPlateEnergy() + partEnergy) / Sim->dynamics.units().unitEnergy()
	  << " " << partpartEnergyLoss  * Sim->dynamics.units().unitTime() / (getTickerTime() * Sim->dynamics.units().unitEnergy())
	  << "\n";
  
  partpartEnergyLoss = 0.0;
}

void 
OPPlateMotion::operator<<(const XMLNode& XML)
{
  try {
    plateName = std::string(XML.getAttribute("PlateName"));
  } catch(...)
    {
      D_throw() << "Could not find the PlateName for the PlateMotion plugin. Did you specify one?";
    }
}

void 
OPPlateMotion::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("PlateMotion");
 
  for (size_t ID(0); ID < localEnergyLoss.size(); ++ID)
    {
      {
	std::fstream of((Sim->dynamics.getLocals()[ID]->getName() 
		       + std::string("EnergyLoss.out")).c_str(), 
		      std::ios::out | std::ios::trunc);
      
	size_t step(0);
	Iflt deltat(getTickerTime() / Sim->dynamics.units().unitTime());
	Iflt sum = localEnergyLoss[ID].first;
	BOOST_FOREACH(const Iflt& val, localEnergyLoss[ID].second)
	  {
	    sum += val;
	    of << deltat * (step++) << " " 
	       << val / (deltat *  Sim->dynamics.units().unitEnergy()) 
	       << "\n";
	  }

	XML << xmlw::tag("Plate")
	    << xmlw::attr("ID") << ID
	    << xmlw::attr("PowerLossRate") 
	    << (sum * Sim->dynamics.units().unitTime() 
		/ (Sim->dSysTime * Sim->dynamics.units().unitEnergy()))
	    << xmlw::endtag("Plate");
      }
      {
	std::fstream of((Sim->dynamics.getLocals()[ID]->getName() 
			 + std::string("EnergyFlux.out")).c_str(), 
			std::ios::out | std::ios::trunc);
	
	size_t step(0);
	Iflt deltat(getTickerTime() / Sim->dynamics.units().unitTime());
	
	BOOST_FOREACH(const Iflt& val, localEnergyFlux[ID].second)
	  of << deltat * (step++) << " " << val / (deltat * Sim->dynamics.units().unitEnergy()) << "\n";
      }
    }
  XML << xmlw::endtag("PlateMotion");
}

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

COPPlateMotion::COPPlateMotion(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPTicker(tmp,"PlateMotion")
{
  operator<<(XML);
}

COPPlateMotion::COPPlateMotion(const COPPlateMotion& cp):
  COPTicker(cp),
  plateID(cp.plateID),
  plateName(cp.plateName)
{
  if (cp.logfile.is_open())
    cp.logfile.close();
}

void
COPPlateMotion::initialise()
{
  try {
    plateID = Sim->Dynamics.getLocal(plateName)->getID();
  } catch(...)
    {
      D_throw() << "Could not find the PlateName specified. You said " << plateName;
    }
  
  if (dynamic_cast<const CLOscillatingPlate*>(Sim->Dynamics.getLocals()[plateID].get_ptr()) == NULL) 
    D_throw() << "The PlateName'd local is not a CLOscillatingPlate";

  if (logfile.is_open())
    logfile.close();
  
  logfile.open("plateMotion.out", std::ios::out|std::ios::trunc);

  localEnergyLoss.resize(Sim->Dynamics.getLocals().size(), std::make_pair(Iflt(0),std::vector<Iflt>()));

  ticker();
}

void 
COPPlateMotion::eventUpdate(const CLocalEvent& localEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    localEnergyLoss[localEvent.getLocalID()].first += pData.getDeltaKE();  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    localEnergyLoss[localEvent.getLocalID()].first += pData.particle1_.getDeltaKE()
    + pData.particle2_.getDeltaKE();
}

void 
COPPlateMotion::ticker()
{
  BOOST_FOREACH(localEntry& entry, localEnergyLoss)
    {
      entry.second.push_back(entry.first);
      entry.first = 0.0;
    }

  Vector com(0,0,0), momentum(0,0,0);
  Iflt sqmom(0);
  
  Iflt mass(0);
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      Vector pos(part.getPosition()), vel(part.getVelocity());
      Iflt pmass(Sim->Dynamics.getSpecies(part).getMass());
      Sim->Dynamics.BCs().applyBC(pos, vel);
      momentum += vel * pmass;
      sqmom += (vel | vel) * (pmass * pmass);
      com += pos * pmass;
      mass += pmass;
    }
  
  com /= (mass * Sim->Dynamics.units().unitLength());
  Vector comvel = momentum / (mass * Sim->Dynamics.units().unitVelocity());
  

  const CLOscillatingPlate& plate(*dynamic_cast<const CLOscillatingPlate*>(Sim->Dynamics.getLocals()[plateID].get_ptr()));

  Vector platePos = (plate.getPosition() - plate.getCentre()) / Sim->Dynamics.units().unitLength();

  Vector plateSpeed = plate.getVelocity() / Sim->Dynamics.units().unitVelocity();

  logfile << Sim->dSysTime / Sim->Dynamics.units().unitTime()
	  << " " << platePos[0] << " " << platePos[1] << " " << platePos[2] 
	  << " " << com[0] << " " << com[1] << " " << com[2]
	  << " " << comvel[0] << " " << comvel[1] << " " << comvel[2]
	  << " " << plateSpeed[0] << " " << plateSpeed[1] << " " << plateSpeed[2]
	  << " " << (sqmom - ((momentum | momentum) / Sim->lN)) / (Sim->lN * pow(Sim->Dynamics.units().unitMomentum(),2))
	  << " " << plate.getPlateEnergy() / Sim->Dynamics.units().unitEnergy()
	  << "\n";
}

void 
COPPlateMotion::operator<<(const XMLNode& XML)
{
  try {
    plateName = std::string(XML.getAttribute("PlateName"));
  } catch(...)
    {
      D_throw() << "Could not find the PlateName for the PlateMotion plugin. Did you specify one?";
    }
}

void 
COPPlateMotion::output(xmlw::XmlStream& XML)
{
  for (size_t ID(0); ID < localEnergyLoss.size(); ++ID)
    {
      std::fstream of((Sim->Dynamics.getLocals()[ID]->getName() 
		       + std::string("EnergyLoss.out")).c_str(), 
		      std::ios::out | std::ios::trunc);
      
      size_t step(0);
      Iflt deltat(getTickerTime() / Sim->Dynamics.units().unitTime());

      BOOST_FOREACH(const Iflt& val, localEnergyLoss[ID].second)
	of << deltat * (step++) << " " << val /Sim->Dynamics.units().unitEnergy() << "\n";
    }
}

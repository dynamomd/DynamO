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

  ticker();
}

void 
COPPlateMotion::ticker()
{
  
  Vector com(0,0,0);
  Iflt mass(0);
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      Vector pos(part.getPosition());
      Iflt pmass(Sim->Dynamics.getSpecies(part).getMass());

      Sim->Dynamics.BCs().setPBC(pos);
      com += pos * pmass;
      mass += pmass;
    }
  
  com /= (mass * Sim->Dynamics.units().unitLength());

  const CLOscillatingPlate& plate(*dynamic_cast<const CLOscillatingPlate*>(Sim->Dynamics.getLocals()[plateID].get_ptr()));

  Vector platePos = (plate.getPosition() - plate.getCentre()) / Sim->Dynamics.units().unitLength();

  logfile << Sim->dSysTime / Sim->Dynamics.units().unitTime()
	  << " " << com[0] << " " << com[1] << " " << com[2]
	  << " " << platePos[0] << " " << platePos[1] << " " << platePos[2] << "\n";
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

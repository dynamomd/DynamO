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

#include "inputplugin.hpp"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/foreach.hpp>
#include "../simulation/particle.hpp"
#include "../schedulers/include.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/species/species.hpp"
#include "../dynamics/interactions/include.hpp"
#include "../dynamics/ranges/include.hpp"
#include "../dynamics/BC/include.hpp"
#include "../dynamics/liouvillean/include.hpp"
#include "../dynamics/systems/ghost.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/topology/include.hpp"

CInputPlugin::CInputPlugin(dynamo::SimData* tmp, const char *aName):
  SimBase(tmp, aName)
{}

void 
CInputPlugin::rescaleVels(double val)
{
  dout << "WARNING Rescaling kT to " << val << std::endl;
  
  double currentkT(Sim->dynamics.getLiouvillean().getkT()
		 / Sim->dynamics.units().unitEnergy());

  dout << "Current kT " << currentkT << std::endl;

  Sim->dynamics.getLiouvillean().rescaleSystemKineticEnergy(val/ currentkT);
}

void 
CInputPlugin::setCOMVelocity(const Vector vel)
{
  dout << "Setting COM Velocity" << std::endl;
  
  if (Sim->N <= 1)
    derr << "Refusing to set momentum for a " 
	     << Sim->N << " particle system" << std::endl;
  else
    Sim->dynamics.setCOMVelocity(vel);  
}


void 
CInputPlugin::zeroMomentum()
{
  dout << "Zeroing Momentum" << std::endl;
  
  if (Sim->N <= 1)
    derr << "Refusing to zero momentum for a " 
	     << Sim->N << " particle system" << std::endl;
  else
    Sim->dynamics.setCOMVelocity();
}

void 
CInputPlugin::zeroCentreOfMass()
{
  dout << "Zeroing Centre of Mass" << std::endl;
  
  Vector com(0,0,0);  
  double totmass = 0.0;
  BOOST_FOREACH(Particle& part, Sim->particleList)  
    {
      totmass += Sim->dynamics.getSpecies(part).getMass(part.getID());
      com += part.getPosition() * Sim->dynamics.getSpecies(part).getMass(part.getID());
    }
  com /= totmass;
  
  BOOST_FOREACH(Particle& part, Sim->particleList)
    part.getPosition() -= com;
}

void 
CInputPlugin::mirrorDirection(unsigned int iDim)
{
  BOOST_FOREACH(Particle& part, Sim->particleList)  
    {
      part.getVelocity()[iDim] *= -1.0;
      part.getPosition()[iDim] *= -1.0;
    }
}

void 
CInputPlugin::zeroVelComp(size_t iDim)
{
  dout << "Zeroing the " << iDim << " dimension velocities" << std::endl;
  BOOST_FOREACH(Particle& part, Sim->particleList)
    part.getVelocity()[iDim] = 0.0;
}

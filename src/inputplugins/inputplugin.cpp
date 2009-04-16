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

#include "inputplugin.hpp"
#include "../extcode/include/boost/random/normal_distribution.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/foreach.hpp>
#include "../simulation/particle.hpp"
#include "../schedulers/include.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/species/species.hpp"
#include "../dynamics/units/include.hpp"
#include "../dynamics/interactions/include.hpp"
#include "../dynamics/ranges/include.hpp"
#include "../dynamics/BC/include.hpp"
#include "../dynamics/liouvillean/include.hpp"
#include "../dynamics/systems/ghost.hpp"
#include "../base/is_exception.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/topology/include.hpp"

CInputPlugin::CInputPlugin(DYNAMO::SimData* tmp, const char *aName, 
			   const char *aColor):
  SimBase(tmp, aName, aColor)
{}

void 
CInputPlugin::rescaleVels(Iflt val)
{
  I_cout() << "WARNING Rescaling kT to " << val;
  
  Iflt currentkT(Sim->Dynamics.Liouvillean().getkT()
		 / Sim->Dynamics.units().unitEnergy());

  I_cout() << "Current kT " << currentkT;

  Vector energy = Sim->Dynamics.Liouvillean().getVectorSystemKineticEnergy();

  Iflt avg  = energy[0];

  for (size_t iDim(1); iDim < NDIM; ++iDim)
    avg += energy[iDim];

  avg /= NDIM;

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    energy[iDim] = sqrt(avg / energy[iDim]);

  Sim->Dynamics.Liouvillean().rescaleSystemKineticEnergy(energy);

  Sim->Dynamics.Liouvillean().rescaleSystemKineticEnergy(val/ currentkT);
}

void 
CInputPlugin::zeroMomentum()
{
  I_cout() << "Zeroing Momentum";    
  Sim->Dynamics.zeroMomentum(Sim->vParticleList);
}

void 
CInputPlugin::zeroCentreOfMass()
{
  I_cout() << "Zeroing Centre of Mass";
  
  Vector com(0,0,0);  
  Iflt totmass = 0.0;
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)  
    {
      totmass += Sim->Dynamics.getSpecies(part).getMass();
      com += part.getPosition() * Sim->Dynamics.getSpecies(part).getMass();
    }
  com /= totmass;
  
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    part.getPosition() -= com;
}

void 
CInputPlugin::setPackFrac(Iflt tmp)
{
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp->getIntPtr()->hardCoreDiam(), NDIM) * sp->getCount();
  
  volume *= PI / (6 * (Sim->Dynamics.units().simVolume()));

  Sim->Dynamics.rescaleLengths(pow(tmp/volume, 1.0/3.0) -1.0);
}

void 
CInputPlugin::mirrorDirection(unsigned int iDim)
{
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)  
    {
      part.getVelocity()[iDim] *= -1.0;
      part.getPosition()[iDim] *= -1.0;
    }
}

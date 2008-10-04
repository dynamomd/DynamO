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
  Iflt energy = 0.0;  
  
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)  
    energy += (part.getVelocity()).square() 
    * Sim->Dynamics.getSpecies(part).getMass();
  
  I_cout() << "WARNING Rescaling kT to " << val;
  Iflt rescale = sqrt (val * 3.0 * Sim->vParticleList.size() * Sim->Dynamics.units().unitEnergy() / energy);
  
  CVector<> vel; //Passes by reference
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    part.getVelocity() = part.getVelocity() * rescale;
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
  
  CVector<> com(0.0);  
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
  
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp.getIntPtr()->hardCoreDiam(), NDIM) * sp.getCount();
  
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

/*
void 
CInputPlugin::setSimType(unsigned int i)
{
  if (!diamScale) D_throw() << "Diamscale not set!";
  
  switch (i)
    {
    case 1:
      //Just a hard sphere system
      Sim->ptrScheduler = new CSMultList(Sim);
      Sim->Dynamics.setPBC<CSPBC>();
      Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
      Sim->Dynamics.addInteraction(new CIHardSphere(Sim, diamScale, 1.0, 
						    new C2RAll()
						    ))->setName("Bulk");
      Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(0,nParts-1), 1.0, "Bulk", 0, "Bulk"));
      Sim->Dynamics.setUnits(new CUElastic(diamScale,Sim));
      break;
    case 2:
      //Just a square well system
      Sim->ptrScheduler = new CSMultList(Sim);
      Sim->Dynamics.setPBC<CSPBC>();
      Sim->Dynamics.setUnits(new CUSW(diamScale,1.0, Sim));
      Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
      Sim->Dynamics.addInteraction(new CISquareWell(Sim, diamScale, 1.5, 1.0, 1.0,
						    new C2RAll()
						    ))->setName("Bulk");
      Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(0,nParts-1), 1.0, "Bulk", 0, "Bulk"));
      break;
    case 3:
      //For the MaGee random walk system
      break;
    case 4:
      //For the binary hard sphere system
      {
	Iflt molfracA = 0.001;
	unsigned long as = 0;
	unsigned long ae = (static_cast<unsigned long>(molfracA * nParts)) - 1;
	unsigned long bs = ae + 1; 
	unsigned long be = nParts - 1;
	
	Sim->Dynamics.setPBC<CSPBC>();
	Sim->ptrScheduler = new CSMultList(Sim);
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim,1.0 * diamScale,1.0,new C2RSingle(new CRRange(as,ae))
			    ))->setName("AA");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim,0.55*diamScale, 1.0,new C2RPair(new CRRange(as,ae),
								new CRRange(bs,be))
			    ))->setName("AB");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim,0.1*diamScale, 1.0,new C2RSingle(new CRRange(bs,be))))->setName("BB");
	
	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(as,ae), 1.0, "AA", 0, "AA"));
	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(bs,be), 1.0, "BB", 0, "BB"));
	Sim->Dynamics.setUnits(new CUElastic(diamScale, Sim));
      }
      break;
    default:
      D_throw() << "Unrecognised sim type";
    }
}


CVector<> 
CInputPlugin::getRandVelVec()
{
  //See http://mathworld.wolfram.com/SpherePointPicking.html
  boost::normal_distribution<Iflt> normdist(0.0, (1.0 / sqrt(NDIM)));
  
  boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution<Iflt> >
    normal_sampler(Sim->ranGenerator, normdist);
  
  CVector<> tmpVec;
  for (int iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_sampler();
  
  return tmpVec;
}

*/

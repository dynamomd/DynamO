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

#include "rings.hpp"
#include <cmath>
#include "../simulation/particle.hpp"
#include "../base/is_simdata.hpp"
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
//#include "../dynamics/globals/ghost.hpp"
#include "../base/is_exception.hpp"
#include "../base/is_simdata.hpp"

CIRings::CIRings(Iflt dens, CVector<long> uC, long nCL, DYNAMO::SimData* tmp):  
  CInputPlugin(tmp, "HomopolymerRings"), 
  density(dens), volume(1.0), cells(uC), 
  maxdim(0),ncells(1),
  chainlength(nCL)
{
  int maxof2 = 0;
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      ncells *= cells[iDim]; 
      if (cells[iDim] > cells[maxdim])
	maxdim = iDim;

      if ((cells[iDim] > cells[maxdim]) && (iDim < 2))
	maxof2 = iDim;
    }
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      aspectRatio[iDim] = ((Iflt) cells[iDim])/((Iflt) cells[maxdim]);
      volume *= aspectRatio[iDim];
    }

  atomDiam = 1.0;
  bondRadius = 0.5;
  bondwidth = 0.025;

  siteAngle = 2.0 * PI / chainlength;
  siteRadius = (bondRadius + 0.25 * bondwidth) / sin(0.5 * siteAngle);
  latticeWidth = 2.0 * (siteRadius + 0.5 * atomDiam) / dens;

  systemWidth = latticeWidth * cells[maxof2];
  
  if ((maxdim > 2) && (cells[maxdim] * atomDiam / dens > systemWidth))
    systemWidth = cells[maxdim] * atomDiam / dens;

  diamScale = 1.0 / systemWidth;
  nParts = ncells * chainlength;

  CVector<> pos(0);
  Iflt rval = siteRadius * diamScale;

  for (unsigned int i = 0; i < chainlength; i++)
    {
      pos[0] = rval * cos(i * siteAngle);
      pos[1] = -rval * sin(i * siteAngle);
      sitelists.push_back(pos);
    }
}

void
CIRings::initialise()
{  
  I_cout() << "No. of particles = " << nParts;
  
  CVector<> position;
  CVector<long> iterVec(0);

  unsigned long nParticles = 0;

  while (iterVec[NDIM - 1] != cells[NDIM - 1])
    {
      for (int iDim = 0; iDim < NDIM; iDim++)
	position[iDim] = ((1.0/cells[iDim]) * ((Iflt) iterVec[iDim]) - 0.5 + (0.5 * (1.0/cells[iDim]) ));
      
      BOOST_FOREACH(CVector<> npos, sitelists)
	Sim->vParticleList.push_back(CParticle(position + npos, getRandVelVec(), nParticles++)); 

      //Now update the displacement vector
      iterVec[0]++;
      
      for (int iDim = 1; iDim < NDIM; iDim++)
        {
          //This increments the next dimension along when
          if (iterVec[iDim - 1] == cells[iDim -1])
            {
              iterVec[iDim - 1] = 0;
              iterVec[iDim] += 1;
            }
        }
    }
  
  rescaleVels();
  zeroMomentum();
}

void 
CIRings::setSimType(unsigned int)
{
  Sim->ptrScheduler = new CSMultList(Sim);
  //Sim->ptrScheduler = new CSFastSingle(Sim);
  Sim->Dynamics.setPBC<CSPBC>();
  Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

  CInteraction* bonds = Sim->Dynamics.addInteraction
    (new CISquareBond(Sim, 2.0 * bondRadius * diamScale,
		      1.0 + (bondwidth / bondRadius), new C2RRings(0,Sim->vParticleList.size()-1, chainlength)));
  
  bonds->setName("Bonds");
  Sim->Dynamics.addInteraction(new CIHardSphere(Sim, diamScale, 1.0, 
						new C2RPair(new CRAll(Sim),
							    new CRAll(Sim))
						))->setName("Bulk");
  
  Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(0,nParts-1), 1.0, "Bulk", 0, "Bulk"));
  Sim->Dynamics.setUnits(new CUElastic(diamScale,Sim));
  
  CVector<> start(0.0), tmp(0.0);

  /*C2RRangeList* range = dynamic_cast<C2RRangeList*>(bonds->getRange().get_ptr());
  
  for (int i = 0; i < ncells; i++)
    {
	for (unsigned int j = 0; j < chainlength - 1; j++)
	range->addPair(i * chainlength + j, i * chainlength + j + 1);
	
	if (chainlength > 2)
	range->addPair(i * chainlength, i * chainlength + chainlength - 1);

      range->addRange(new C2RRing(i*chainlength, i * chainlength + chainlength - 1));
    }*/

  rescaleVels();
  zeroMomentum();

}

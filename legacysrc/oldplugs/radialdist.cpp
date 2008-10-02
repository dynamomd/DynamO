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

#include "radialdist.hpp"
#include "../dynamics/include.hpp"
#include "../base/is_simdata.hpp"

COPRadialDist::COPRadialDist(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"RadialDist"),
  sampleCount(0)
{
  binWidth = Sim->Dynamics.units().unitLength() * widthScale;

  //+1 for rounding truncation, +1 for the zero bin
  binCount = 2 + (long) (0.5 / binWidth);

  bin.resize(binCount, 0);
}

COPRadialDist::~COPRadialDist()
{}

void
COPRadialDist::collisionUpdate(const CIntEvent &, const CIntEventData &preColl)
{
  sampleCount++;
  CVector<> centerParticle, rij;
  
  for (std::vector<CParticle>::const_iterator iPtr = Sim->vParticleList.begin();
       iPtr != Sim->vParticleList.end();)
    {
      centerParticle = iPtr->getPosition();
      iPtr++;
      
      for (std::vector<CParticle>::const_iterator jPtr = iPtr;
	   jPtr != Sim->vParticleList.end(); jPtr++)
	{
	  rij = centerParticle - jPtr->getPosition();
	  Sim->Dynamics.BCs().setPBC(rij);
	  
	  int i = (long) (((rij.length())/binWidth) + 0.5);

	  //Some particles are further than 0.5 away
	  if (i < binCount)
	    bin[i]++;
	}
    }
}


void
COPRadialDist::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("RadialDistribution") << xmlw::chardata();
  
  long originsTaken = (sampleCount * Sim->vParticleList.size())/2,
    maxinnershells = (long) (1 / widthScale) + 1;

  Iflt radius, volshell, GR, density = getNumberDensity();

  for (long i = 1; i < binCount; i++)
    {
      if (bin[i] > 0)
      {
	radius = widthScale * i;
	volshell = (4.0 * PI * widthScale * radius * radius) +
	((PI * widthScale * widthScale * widthScale) / 3.0);
	GR = ((Iflt) bin[i]) / (density * originsTaken * volshell);
	if (i < maxinnershells)
	{
	  radius = 1.0 + 0.5 * widthScale;
	  GR *= 2.0;
	}
	XML << " " << radius << " " << GR << "\n"; 
      }
    }
  
  XML << xmlw::endtag("RadialDistribution");
}

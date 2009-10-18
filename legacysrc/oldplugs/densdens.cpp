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

#include "densdens.hpp"
#include <boost/foreach.hpp>
#include "../datatypes/complex.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../dynamics/include.hpp"
#include "../base/is_simdata.hpp"

OPDens::OPDens(DYNAMO::SimData* tmp):
  OutputPlugin(tmp, "DensityCorrelation"),
  counter(0),
  collcount(0)
{}

OPDens::~OPDens()
{}

void 
OPDens::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  collcount++;

  if ((collcount % 1000) == 0)
    for (std::vector<CParticle>::const_iterator iPtr = Sim->vParticleList.begin();
	 iPtr != Sim->vParticleList.end(); iPtr++)
      {
	counter++;
	
	Vector  pos = iPtr->getPosition();
	
	for (long qy = 0; qy < 20; qy++)
	  for (long qz = 0; qz < 20; qz++)
	    xybin[qy][qz] += CComplex(0.0, 2.0*PI*(((Iflt) qy)*pos[1]+((Iflt)qz)*pos[2])).exponent();
	
	/*for (std::vector<CParticle>::const_iterator jPtr = iPtr + 1;
	  jPtr != sim.particleList.end(); jPtr++)
	  {
	  pos = iPtr->getPosition() - jPtr->getPosition();
	  for (long qy = 0; qy < 20; qy++)
	  for (long qz = 0; qz < 20; qz++)
	  densdens[qy][qz] += CComplex(0.0, 2.0*PI*(((Iflt)qy)*pos[1]+((Iflt)qz)*pos[2])).exponent();
	  }*/
      }
  
}

void
OPDens::output(xmlw::XmlStream &XML)
{
  Iflt factor = 1.0 / ((Iflt) Sim->vParticleList.size() * counter);
  
  XML << xmlw::tag("rho(q)") << xmlw::chardata();
  for (long qy = 0; qy < 20; qy++)
    for (long qz = 0; qz < 20; qz++)
      XML << qy << " " << qz << " " << xybin[qy][qz].getr()*factor 
	  << " " << xybin[qy][qz].geti()*factor << "\n";
  XML << xmlw::endtag("rho(q)")
      << xmlw::tag("|rho(q)|^2") << xmlw::chardata();

  for (long qy = 0; qy < 20; qy++)
    for (long qz = 0; qz < 20; qz++)
      XML << qy << " " << qz << " " 
	  << pow(xybin[qy][qz].geti()*factor, 2.0) 
	+ pow(xybin[qy][qz].getr()*factor, 2.0) 
	  << "\n";

  XML << xmlw::endtag("|rho(q)|^2");
}

void
OPDens::periodicOutput()
{}

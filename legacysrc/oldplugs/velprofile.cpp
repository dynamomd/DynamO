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

#include "velprofile.hpp"
#include <boost/foreach.hpp>
#include "../extcode/xmlwriter.hpp"
#include "../datatypes/vector.hpp"
#include "../dynamics/include.hpp"
#include "../base/is_simdata.hpp"

#define NBins 20

COPVPROF::COPVPROF(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"VelProfile"),
  vxy(1.0/NBins, -0.5, NBins),
  counter(1.0/NBins, -0.5, NBins),
  samplesTaken(0)
{}

COPVPROF::~COPVPROF()
{}

void 
COPVPROF::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  CVector<> pos;
  samplesTaken++;

  if (!(Sim->lNColl % (Sim->vParticleList.size()/10)))
    BOOST_FOREACH( const CParticle & Part, Sim->vParticleList )
    {
      pos = Part.getPosition();
      Sim->Dynamics.BCs().setPBC(pos);
      vxy[pos[1]][pos[2]] += (Sim->Dynamics.getLabVelocity(Part))[0];
      counter[pos[1]][pos[2]]++;
    }
}

void
COPVPROF::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Vx_y_velprofile")
      << xmlw::chardata();

  for (long y = 0; y < NBins; y++)
    for (long z = 0; z < NBins; z++)
      if (counter[y][z] != 0)
	XML << y*counter.binWidth - 0.5 << " " << z*counter.binWidth - 0.5 
	    << " " << vxy[y][z]/((Iflt) counter[y][z]) << "\n";
      else
	XML << y << " " << z << " 0\n"; 
  
  XML << xmlw::endtag("Vx_y_velprofile")
      << xmlw::tag("xy_density")
      << xmlw::chardata();

  for (long y = 0; y < NBins; y++)
    for (long z = 0; z < NBins; z++)
      if (counter[y][z] != 0)
	XML << y*counter.binWidth - 0.5 << " " << z*counter.binWidth - 0.5 
	    << " " << ((Iflt) counter[y][z]) << "\n";
      else
	XML << y << " " << z << " 0\n"; 
  
  XML << xmlw::endtag("xy_density");

}

void
COPVPROF::periodicOutput()
{}

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

#include "mft.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../dynamics/include.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../base/is_simdata.hpp"

COPMFT::COPMFT(DYNAMO::SimData* tmp):
  COutputPlugin(tmp, "MeanFreeTime"),
  mft(0.01*Sim->Dynamics.units().unitLength())
{
  //Initialise the array
  for (size_t i = 0; i < Sim->vParticleList.size(); i++)
    particle2time[Sim->vParticleList[i].getID()] = 0.0;
}

COPMFT::~COPMFT()
{}

void 
COPMFT::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  if (particle2time[collision.getParticle1().getID()] != 0.0)
    mft.addVal(Sim->dSysTime - particle2time[collision.getParticle1().getID()]);

  if (particle2time[collision.getParticle2().getID()] != 0.0)
    mft.addVal(Sim->dSysTime - particle2time[collision.getParticle2().getID()]);

  particle2time[collision.getParticle1().getID()] = Sim->dSysTime;
  particle2time[collision.getParticle2().getID()] = Sim->dSysTime;
}

void
COPMFT::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("MFT");
  mft.outputHistogram(XML, 1.0/Sim->Dynamics.units().unitTime());
  XML << xmlw::endtag("MFT");
}

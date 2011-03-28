/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "collisionCorrelator.hpp"
#include "../../dynamics/2particleEventData.hpp"
#include "../../dynamics/include.hpp"
#include <magnet/xmlwriter.hpp>

OPCollisionCorrelator::OPCollisionCorrelator(const DYNAMO::SimData* t1,
					       const XMLNode& XML):
  OP2PP(t1,"CollisionCorrelator")
{ operator<<(XML); }

void 
OPCollisionCorrelator::operator<<(const XMLNode& XML)
{
  try {
    /*if (XML.isAttributeSet("length"))
      collisionHistoryLength 
	= boost::lexical_cast<size_t>
	(XML.getAttribute("length"));*/
      }
  catch (std::exception& excep)
    {
      M_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void 
OPCollisionCorrelator::initialise()
{
  //Set the history size
  lastColl.resize(Sim->N, std::vector<double>(Sim->N, 0.0));

  if (Sim->lastRunMFT == 0.0)
    M_throw() << "This output plugin requires an estimate for the mean free time. run the configuration a little first.";

  //Histogram in mean free times
  freetimehist = C1DHistogram(Sim->lastRunMFT*0.1);
}

void 
OPCollisionCorrelator::A2ParticleChange(const PairEventData& PDat)
{
  size_t ID1 = PDat.particle1_.getParticle().getID();
  size_t ID2 = PDat.particle2_.getParticle().getID();

  if (ID1 > ID2) std::swap(ID1, ID2);

  //Check there was a previous collision
  if (lastColl[ID1][ID2] != 0.0)
    freetimehist.addVal(Sim->dSysTime - lastColl[ID1][ID2]);

  lastColl[ID1][ID2] = Sim->dSysTime;  
}

void
OPCollisionCorrelator::output(xml::XmlStream &XML)
{
  for (size_t ID1(0); ID1 < Sim->N; ++ID1)
    for (size_t ID2(ID1+1); ID2 < Sim->N; ++ID2)
      if (lastColl[ID1][ID2] > 100* freetimehist.data.binWidth) freetimehist.addVal(-1.0);

  XML << xml::tag("CollisionCorrelator");
  
  freetimehist.outputHistogram(XML, 1.0/Sim->dynamics.units().unitTime());

  XML << xml::endtag("CollisionCorrelator");
}

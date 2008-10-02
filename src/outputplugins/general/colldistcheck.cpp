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

#include "colldistcheck.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../0partproperty/collMatrix.hpp"

COPCollDistCheck::COPCollDistCheck(const DYNAMO::SimData* t1):
  COutputPlugin(t1,"CollDistCheck"),
  ptrCM(NULL)
{
}

COPCollDistCheck::~COPCollDistCheck()
{
}

void 
COPCollDistCheck::initialise() 
{
  ptrCM = Sim->getOutputPlugin<COPCollMatrix>();
}

void 
COPCollDistCheck::eventUpdate(const CIntEvent& eevent, 
			      const C2ParticleData& PDat)
{
  std::pair<size_t, EEventType> locPair(ptrCM->getID(eevent.getInteraction()), 
					eevent.getType());

  if (distList.find(locPair) == distList.end())
    distList[locPair] = C1DHistogram(0.01 * Sim->Dynamics.units().unitLength());
  
  distList[locPair].addVal(PDat.rij.length());
}

void 
COPCollDistCheck::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("CollDistCheck");
  
  typedef std::pair<const std::pair<size_t, EEventType>, C1DHistogram>
    localPair;

  BOOST_FOREACH(const localPair& p, distList)
    {
      XML << xmlw::tag("Distance") << xmlw::attr("Name") 
	  << ptrCM->getName(p.first.first)
	  << xmlw::attr("EventType") 
	  << CIntEvent::getCollEnumName(p.first.second);

      p.second.outputHistogram(XML, 1.0 / Sim->Dynamics.units().unitLength());

      
      XML << xmlw::endtag("Distance");
    }
  
  XML << xmlw::endtag("CollDistCheck");

}

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

#include "colldistcheck.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"

OPCollDistCheck::OPCollDistCheck(const DYNAMO::SimData* t1, 
				   const XMLNode& XML):
  OutputPlugin(t1,"CollDistCheck"),
  binwidth(0.01)
{
  operator<<(XML);
}

OPCollDistCheck::~OPCollDistCheck()
{
}

void 
OPCollDistCheck::initialise() 
{}

void 
OPCollDistCheck::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth"))
	binwidth = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPCorrelator";
    }
  
}

void 
OPCollDistCheck::eventUpdate(const CIntEvent& eevent, 
			      const C2ParticleData& PDat)
{
  const eventKey locPair(getClassKey(eevent), 
			 eevent.getType());
  
  if (distList.find(locPair) == distList.end())
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
 
  distList[locPair].addVal(PDat.rij.nrm());
}

void 
OPCollDistCheck::eventUpdate(const CGlobEvent& gEvent, 
			      const CNParticleData& PDat)
{
  const eventKey locPair(getClassKey(gEvent), 
			 gEvent.getType());
  
  if ((!PDat.L2partChanges.empty()) 
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const C2ParticleData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}

void 
OPCollDistCheck::eventUpdate(const CLocalEvent& lEvent, 
			      const CNParticleData& PDat)
{
  const eventKey locPair(getClassKey(lEvent), 
			 lEvent.getType());
  
  if ((!PDat.L2partChanges.empty()) 
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const C2ParticleData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}
  
void 
OPCollDistCheck::eventUpdate(const CSystem& sysEvent, 
			      const CNParticleData& PDat, const Iflt& dt)
{
  const eventKey locPair(getClassKey(sysEvent), sysEvent.getType());
  
  if ((!PDat.L2partChanges.empty())
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const C2ParticleData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}

void 
OPCollDistCheck::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("CollDistCheck");
  
  typedef std::pair<eventKey, C1DHistogram>
    localPair;
  
  BOOST_FOREACH(const localPair& p, distList)
    {
      XML << xmlw::tag("Distance") << xmlw::attr("Name") 
	  << getName(p.first.first, Sim)
	  << xmlw::attr("Type") 
	  << CIntEvent::getCollEnumName(p.first.first.second)
	  << xmlw::attr("EventType") 
	  << CIntEvent::getCollEnumName(p.first.second);

      p.second.outputHistogram(XML, 1.0 / Sim->dynamics.units().unitLength());
      
      XML << xmlw::endtag("Distance");
    }
  
  XML << xmlw::endtag("CollDistCheck");

}

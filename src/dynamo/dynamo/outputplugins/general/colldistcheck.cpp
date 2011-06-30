/*  dynamo:- Event driven molecular dynamics simulator 
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

#include "colldistcheck.hpp"
#include "../../dynamics/include.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

OPCollDistCheck::OPCollDistCheck(const dynamo::SimData* t1, 
				 const magnet::xml::Node& XML):
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
OPCollDistCheck::operator<<(const magnet::xml::Node& XML)
{
  try 
    {
      binwidth = XML.getAttribute("binwidth").as<double>(0.01);
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPCorrelator";
    }
  
}

void 
OPCollDistCheck::eventUpdate(const IntEvent& eevent, 
			      const PairEventData& PDat)
{
  const eventKey locPair(getClassKey(eevent), 
			 eevent.getType());
  
  if (distList.find(locPair) == distList.end())
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
 
  distList[locPair].addVal(PDat.rij.nrm());
}

void 
OPCollDistCheck::eventUpdate(const GlobalEvent& gEvent, 
			      const NEventData& PDat)
{
  const eventKey locPair(getClassKey(gEvent), 
			 gEvent.getType());
  
  if ((!PDat.L2partChanges.empty()) 
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const PairEventData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}

void 
OPCollDistCheck::eventUpdate(const LocalEvent& lEvent, 
			      const NEventData& PDat)
{
  const eventKey locPair(getClassKey(lEvent), 
			 lEvent.getType());
  
  if ((!PDat.L2partChanges.empty()) 
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const PairEventData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}
  
void 
OPCollDistCheck::eventUpdate(const System& sysEvent, 
			      const NEventData& PDat, const double& dt)
{
  const eventKey locPair(getClassKey(sysEvent), sysEvent.getType());
  
  if ((!PDat.L2partChanges.empty())
      && (distList.find(locPair) == distList.end()))
    distList[locPair] = C1DHistogram(binwidth * Sim->dynamics.units().unitLength());
  
  BOOST_FOREACH(const PairEventData& dat, PDat.L2partChanges)
    distList[locPair].addVal(dat.rij.nrm());
}

void 
OPCollDistCheck::output(magnet::xml::XmlStream& XML)
{
  XML << magnet::xml::tag("CollDistCheck");
  
  typedef std::pair<eventKey, C1DHistogram>
    localPair;
  
  BOOST_FOREACH(const localPair& p, distList)
    {
      XML << magnet::xml::tag("Distance") << magnet::xml::attr("Name") 
	  << getName(p.first.first, Sim)
	  << magnet::xml::attr("Type") 
	  << p.first.first.second
	  << magnet::xml::attr("EventType") 
	  << p.first.second;

      p.second.outputHistogram(XML, 1.0 / Sim->dynamics.units().unitLength());
      
      XML << magnet::xml::endtag("Distance");
    }
  
  XML << magnet::xml::endtag("CollDistCheck");

}

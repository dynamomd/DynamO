/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/collMatrix.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/include.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPCollMatrix::OPCollMatrix(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"CollisionMatrix"),
    totalCount(0)
  {
  }

  void 
  OPCollMatrix::initialise()
  {
    lastEvent.resize(Sim->N(), lastEventData(Sim->systemTime, eventKey(classKey(0, NOSOURCE), NONE)));
  }

  OPCollMatrix::~OPCollMatrix()
  {}

  void 
  OPCollMatrix::eventUpdate(const Event& event, const NEventData& SDat)
  {
    for (const ParticleEventData& pData : SDat.L1partChanges)
      newEvent(pData.getParticleID(), pData.getType(), getClassKey(event));  
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	newEvent(pData.particle1_.getParticleID(), pData.getType(), getClassKey(event));
	newEvent(pData.particle2_.getParticleID(), pData.getType(), getClassKey(event));
      }
  }

  void 
  OPCollMatrix::newEvent(const size_t& part, const EEventType& etype, const classKey& ck)
  {
    if (lastEvent[part].second.first.second != NOSOURCE)
      {
	counterData& refCount = counters[counterKey(eventKey(ck,etype), lastEvent[part].second)];
      
	refCount.totalTime += Sim->systemTime - lastEvent[part].first;
	++(refCount.count);
	++(totalCount);
      }
    else
      ++initialCounter[eventKey(ck,etype)];

    lastEvent[part].first = Sim->systemTime;
    lastEvent[part].second = eventKey(ck, etype);
  }

  void
  OPCollMatrix::output(magnet::xml::XmlStream &XML)
  {
  
    XML << magnet::xml::tag("CollCounters") 
	<< magnet::xml::tag("TransitionMatrix");
  
    std::map<eventKey, std::pair<size_t, double> > totmap;
  
    typedef std::pair<const counterKey, counterData> locPair;
  
  
    size_t initialsum(0);
  
    for (const auto& n : initialCounter)
      initialsum += n.second;
  
    for (const locPair& ele : counters)
      {
	XML << magnet::xml::tag("Count")
	    << magnet::xml::attr("Event") << ele.first.first.second
	    << magnet::xml::attr("Name") << getName(ele.first.first.first, Sim)
	    << magnet::xml::attr("lastEvent") << ele.first.second.second
	    << magnet::xml::attr("lastName") << getName(ele.first.second.first, Sim)
	    << magnet::xml::attr("Percent") << 100.0 * ((double) ele.second.count) 
	  / ((double) totalCount)
	    << magnet::xml::attr("mft") << ele.second.totalTime
	  / (Sim->units.unitTime() * ((double) ele.second.count))
	    << magnet::xml::endtag("Count");
      
	//Add the total count
	totmap[ele.first.first].first += ele.second.count;
      
	//Add the rate
	totmap[ele.first.first].second += ((double) ele.second.count) 
	  / ele.second.totalTime;
      }
  
    XML << magnet::xml::endtag("TransitionMatrix")
	<< magnet::xml::tag("Totals");
  
    for (const auto& mp1 : totmap)
      XML << magnet::xml::tag("TotCount")
	  << magnet::xml::attr("Name") << getName(mp1.first.first, Sim)
	  << magnet::xml::attr("Event") << mp1.first.second
	  << magnet::xml::attr("Percent") 
	  << 100.0 * (((double) mp1.second.first)
		      +((double) initialCounter[mp1.first]))
      / (((double) totalCount) + ((double) initialsum))
	  << magnet::xml::attr("Count") << mp1.second.first + initialCounter[mp1.first]
	  << magnet::xml::attr("EventMeanFreeTime")
	  << Sim->systemTime / ((mp1.second.first + initialCounter[mp1.first])
			      * Sim->units.unitTime())
	  << magnet::xml::endtag("TotCount");
  
    XML << magnet::xml::endtag("Totals")
	<< magnet::xml::endtag("CollCounters");
  }
}

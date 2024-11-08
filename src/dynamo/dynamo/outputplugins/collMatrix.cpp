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

    for (const auto& p1: Sim->particles)
      for (const auto& p2: Sim->particles)
        if (p1 != p2)
        {
          auto iPtr = std::dynamic_pointer_cast<ICapture>(Sim->getInteraction(p1, p2));
          if (iPtr)
            {
              auto status = iPtr->isCaptured(p1, p2);
              _lastCaptureState[captureKey(iPtr->getID(), p1)] += status;
              _lastCaptureState[captureKey(iPtr->getID(), p2)] += status;
            }
        }  
  }

  OPCollMatrix::~OPCollMatrix()
  {}

  void 
  OPCollMatrix::eventUpdate(const Event& event, const NEventData& SDat)
  {
    auto ck = getClassKey(event);

    for (const ParticleEventData& pData : SDat.L1partChanges)
      newEvent(pData.getParticleID(), pData.getType(), ck);  
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	      newEvent(pData.particle1_.getParticleID(), pData.getType(), ck);
	      newEvent(pData.particle2_.getParticleID(), pData.getType(), ck);
      }

    //Check if the interaction is a subclass of ICapture, if so grab a shared_ptr
    //to the interaction and call the capture function.
    if (event._source == INTERACTION)
      {
        auto iPtr = std::dynamic_pointer_cast<ICapture>(Sim->interactions[event._sourceID]);
        if (iPtr) {
          for (const PairEventData& pData : SDat.L2partChanges) {
            auto ck1 = captureKey(iPtr->getID(), pData.particle1_.getParticleID());
            auto ck2 = captureKey(iPtr->getID(), pData.particle2_.getParticleID());
            auto& cs1 = _lastCaptureState[ck1];
            auto& cs2 = _lastCaptureState[ck2];

            auto ek = eventKey(ck, pData.getType());
            auto cek1 = captureEventKey(ek, cs1);
            auto cek2 = captureEventKey(ek, cs2);
            // Update the tracked capture status 
            _captureCounters[cek1].count += 1;
            _captureCounters[cek2].count += 1;

            // Histogram of capture state occupancy times
            // Time weighted average of capture state

            auto status = iPtr->isCaptured(pData.particle1_.getParticleID(), pData.particle2_.getParticleID());
            // Update the tracked capture status
            if ((pData.getType() == STEP_OUT) && (!status))
              {
                _lastCaptureState[ck1] -= 1;
                _lastCaptureState[ck2] -= 1;
              }
            if ((pData.getType() == STEP_IN) && (status == 1))
              {
                _lastCaptureState[ck1] += 1;
                _lastCaptureState[ck2] += 1;
              }
          }
        }
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
  
    XML << magnet::xml::endtag("Totals");

    XML	<< magnet::xml::tag("CaptureCounters");
    for (const auto& val: _captureCounters) {
      auto cek =  val.first;
      auto ek = cek.first;
      auto class_key = ek.first;
      auto event_type = ek.second;
      auto captures = cek.second;
      auto captureData = val.second;

      XML << magnet::xml::tag("Count")
          << magnet::xml::attr("Name") << getName(class_key, Sim)
          << magnet::xml::attr("Event") << event_type
          << magnet::xml::attr("captures") << captures
          << magnet::xml::attr("Count") << captureData.count
          << magnet::xml::endtag("Count")
      ; 
    }
    XML	<< magnet::xml::endtag("CaptureCounters");

    XML << magnet::xml::endtag("CollCounters");
  }
}

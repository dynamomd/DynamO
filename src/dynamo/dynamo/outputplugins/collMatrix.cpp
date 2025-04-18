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

#include <dynamo/include.hpp>
#include <dynamo/interactions/captures.hpp>
#include <dynamo/outputplugins/collMatrix.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPCollMatrix::OPCollMatrix(const dynamo::Simulation *tmp,
                           const magnet::xml::Node &)
    : OutputPlugin(tmp, "CollisionMatrix"), totalCount(0),
      _captureStateHistogram(1.0) {}

void OPCollMatrix::initialise() {
  lastEvent.clear();
  lastEvent.resize(Sim->N(),
                   lastEventData(Sim->systemTime,
                                 EventKey(EventSourceKey(0, NOSOURCE), NONE)));

  // Reset the current capture state
  _currentCaptureState.clear();
  for (size_t i(0); i < Sim->N(); ++i)
    for (size_t j(i + 1); j < Sim->N(); ++j) {
      auto &p1 = Sim->particles[i];
      auto &p2 = Sim->particles[j];
      auto iPtr =
          std::dynamic_pointer_cast<ICapture>(Sim->getInteraction(p1, p2));
      if (iPtr) {
        auto status = iPtr->isCaptured(p1, p2) > 0;
        _currentCaptureState[TotalCaptureStateKey(iPtr->getID(), p1)]._state +=
            status;
        _currentCaptureState[TotalCaptureStateKey(iPtr->getID(), p2)]._state +=
            status;
      }
    }
  // Move the time origin to now
  for (auto &p : _currentCaptureState)
    p.second._last_update = Sim->systemTime;
}

OPCollMatrix::~OPCollMatrix() {}

void OPCollMatrix::eventUpdate(const Event &event, const NEventData &SDat) {
  auto ck = getEventSourceKey(event);

  // Here we're investigating particles by capture state
  // We do it first, before last_event is updated.
  //
  // Check if the interaction is a subclass of ICapture, if so grab a shared_ptr
  // to the interaction and call the capture function.
  if (event._source == INTERACTION) {
    auto iPtr =
        std::dynamic_pointer_cast<ICapture>(Sim->interactions[event._sourceID]);
    if (iPtr) {
      for (const PairEventData &pData : SDat.L2partChanges) {
        auto ck1 = TotalCaptureStateKey(iPtr->getID(),
                                        pData.particle1_.getParticleID());
        auto ck2 = TotalCaptureStateKey(iPtr->getID(),
                                        pData.particle2_.getParticleID());

        auto &cs1 = _currentCaptureState[ck1];
        auto &cs2 = _currentCaptureState[ck2];

        auto ek = EventKey(ck, pData.getType());
        auto cek1 = EventCaptureStateKey(ek, cs1._state);
        auto cek2 = EventCaptureStateKey(ek, cs2._state);

        // Lookup capture counter data, but initialise histogram if it doesn't
        // exist
        auto it =
            _captureCounters.insert(decltype(_captureCounters)::value_type(
                cek1, EventCaptureStateData(Sim->lastRunMFT * 0.01)));
        auto &cekd1 = it.first->second;
        it = _captureCounters.insert(decltype(_captureCounters)::value_type(
            cek2, EventCaptureStateData(Sim->lastRunMFT * 0.01)));
        auto &cekd2 = it.first->second;

        // We only track the time between events of the same type, at the
        // start of the simulation we don't have a previous event so we skip.
        if (cekd1.last_event_time != 0)
          cekd1.MFT.addVal(Sim->systemTime - cekd1.last_event_time);
        if (cekd2.last_event_time != 0)
          cekd2.MFT.addVal(Sim->systemTime - cekd2.last_event_time);

        cekd1.rijdotvij.addVal(pData.rvdot);
        cekd2.rijdotvij.addVal(pData.rvdot);

        auto rijdotdP = pData.rij * pData.impulse;
        cekd1.rijdotdP.addVal(rijdotdP);
        cekd2.rijdotdP.addVal(rijdotdP);

        cekd1.vi2.addVal(pData.particle1_.getOldVel().nrm2());
        cekd2.vi2.addVal(pData.particle2_.getOldVel().nrm2());
        // Now we update the last event time
        cekd1.last_event_time = Sim->systemTime;
        cekd2.last_event_time = Sim->systemTime;

        auto new_capture_state = iPtr->isCaptured(
            pData.particle1_.getParticleID(), pData.particle2_.getParticleID());

        if (cs1._last_event_time != 0) {
          cekd1._particle_MFT.addVal(Sim->systemTime - cs1._last_event_time);
          MFTKey k1(EventCaptureStateKey(cs1._last_event, cs1._state), cek1);
          auto it2 = _fullMFT.insert(decltype(_fullMFT)::value_type(
              k1, magnet::math::Histogram<>(Sim->lastRunMFT * 0.01)));
          it2.first->second.addVal(Sim->systemTime - cs1._last_event_time);
        }
        if (cs2._last_event_time != 0) {
          cekd2._particle_MFT.addVal(Sim->systemTime - cs2._last_event_time);
          MFTKey k2(EventCaptureStateKey(cs2._last_event, cs2._state), cek2);
          auto it2 = _fullMFT.insert(decltype(_fullMFT)::value_type(
              k2, magnet::math::Histogram<>(Sim->lastRunMFT * 0.01)));
          it2.first->second.addVal(Sim->systemTime - cs2._last_event_time);
        }

        _captureStateHistogram.addVal(cs1._state,
                                      Sim->systemTime - cs1._last_update);
        _captureStateHistogram.addVal(cs2._state,
                                      Sim->systemTime - cs2._last_update);
        cs1._last_update = Sim->systemTime;
        cs2._last_update = Sim->systemTime;

        // Update the tracked capture/event status
        cs1._last_event = ek;
        cs2._last_event = ek;
        cs1._last_event_time = Sim->systemTime;
        cs2._last_event_time = Sim->systemTime;
        if ((pData.getType() == STEP_OUT) && (new_capture_state == 0)) {
          cs1._state -= 1;
          cs2._state -= 1;
        }
        if ((pData.getType() == STEP_IN) && (new_capture_state == 1)) {
          cs1._state += 1;
          cs2._state += 1;
        }
      }
    }
  }

  for (const ParticleEventData &pData : SDat.L1partChanges)
    newEvent(pData.getParticleID(), pData.getType(), ck);

  for (const PairEventData &pData : SDat.L2partChanges) {
    newEvent(pData.particle1_.getParticleID(), pData.getType(), ck);
    newEvent(pData.particle2_.getParticleID(), pData.getType(), ck);
  }
}
void OPCollMatrix::newEvent(const size_t &part, const EEventType &etype,
                            const EventSourceKey &ck) {
  if (lastEvent[part].second.first.second != NOSOURCE) {
    InterEventData &refCount =
        counters[InterEventKey(EventKey(ck, etype), lastEvent[part].second)];

    refCount.totalTime += Sim->systemTime - lastEvent[part].first;
    ++(refCount.count);
    ++(totalCount);
  } else
    ++initialCounter[EventKey(ck, etype)];

  lastEvent[part].first = Sim->systemTime;
  lastEvent[part].second = EventKey(ck, etype);
}

void OPCollMatrix::output(magnet::xml::XmlStream &XML) {

  XML << magnet::xml::tag("CollCounters")
      << magnet::xml::tag("TransitionMatrix");

  std::map<EventKey, std::pair<size_t, double>> totmap;

  typedef std::pair<const InterEventKey, InterEventData> locPair;

  size_t initialsum(0);

  for (const auto &n : initialCounter)
    initialsum += n.second;

  for (const locPair &ele : counters) {
    XML << magnet::xml::tag("Count") << magnet::xml::attr("Event")
        << ele.first.first.second << magnet::xml::attr("Name")
        << getEventSourceName(ele.first.first.first, Sim)
        << magnet::xml::attr("lastEvent") << ele.first.second.second
        << magnet::xml::attr("lastName")
        << getEventSourceName(ele.first.second.first, Sim)
        << magnet::xml::attr("Percent")
        << 100.0 * ((double)ele.second.count) / ((double)totalCount)
        << magnet::xml::attr("mft")
        << ele.second.totalTime /
               (Sim->units.unitTime() * ((double)ele.second.count))
        << magnet::xml::endtag("Count");

    // Add the total count
    totmap[ele.first.first].first += ele.second.count;

    // Add the rate
    totmap[ele.first.first].second +=
        ((double)ele.second.count) / ele.second.totalTime;
  }

  XML << magnet::xml::endtag("TransitionMatrix") << magnet::xml::tag("Totals");

  for (const auto &mp1 : totmap)
    XML << magnet::xml::tag("TotCount") << magnet::xml::attr("Name")
        << getEventSourceName(mp1.first.first, Sim)
        << magnet::xml::attr("Event") << mp1.first.second
        << magnet::xml::attr("Percent")
        << 100.0 *
               (((double)mp1.second.first) +
                ((double)initialCounter[mp1.first])) /
               (((double)totalCount) + ((double)initialsum))
        << magnet::xml::attr("Count")
        << mp1.second.first + initialCounter[mp1.first]
        << magnet::xml::attr("EventMeanFreeTime")
        << Sim->systemTime / ((mp1.second.first + initialCounter[mp1.first]) *
                              Sim->units.unitTime())
        << magnet::xml::endtag("TotCount");

  XML << magnet::xml::endtag("Totals");

  XML << magnet::xml::tag("CaptureCounters");
  for (const auto &val : _captureCounters) {
    auto cek = val.first;
    auto ek = cek.first;
    auto class_key = ek.first;
    auto event_type = ek.second;
    auto captures = cek.second;
    auto EventCaptureStateData = val.second;

    XML << magnet::xml::tag("Count") << magnet::xml::attr("Name")
        << getEventSourceName(class_key, Sim) << magnet::xml::attr("Event")
        << event_type << magnet::xml::attr("captures") << captures;

    XML << magnet::xml::tag("MFT");
    val.second.MFT.outputHistogram(XML, 1.0 / Sim->units.unitTime());
    XML << magnet::xml::endtag("MFT");

    XML << magnet::xml::tag("ParticleMFT");
    val.second._particle_MFT.outputHistogram(XML, 1.0 / Sim->units.unitTime());
    XML << magnet::xml::endtag("ParticleMFT");

    XML << magnet::xml::tag("RijDotVij");
    val.second.rijdotvij.outputHistogram(XML, 1.0 / Sim->units.unitLength() /
                                                  Sim->units.unitVelocity());
    XML << magnet::xml::endtag("RijDotVij");

    XML << magnet::xml::tag("RijDotDeltaPij");
    val.second.rijdotvij.outputHistogram(XML, 1.0 / Sim->units.unitLength() /
                                                  Sim->units.unitMomentum());
    XML << magnet::xml::endtag("RijDotDeltaPij");

    XML << magnet::xml::tag("V2");
    val.second.vi2.outputHistogram(XML, 1.0 / Sim->units.unitVelocity() /
                                            Sim->units.unitVelocity());
    XML << magnet::xml::endtag("V2");

    XML << magnet::xml::endtag("Count");
  }
  XML << magnet::xml::endtag("CaptureCounters")
      << magnet::xml::tag("CaptureStateHistogram");

  // Before we output the histogram we need to bring everything up to date
  for (auto &p : _currentCaptureState) {
    auto &cs = p.second;
    _captureStateHistogram.addVal(cs._state, Sim->systemTime - cs._last_update);
    cs._last_update = Sim->systemTime;
  }

  _captureStateHistogram.outputHistogram(XML, 1.0 / Sim->units.unitEnergy());
  XML << magnet::xml::endtag("CaptureStateHistogram")
      << magnet::xml::endtag("CollCounters");

  XML << magnet::xml::tag("FullMFTs");

  for (auto &p : _fullMFT) {
    auto &k = p.first;
    auto &h = p.second;
    XML << magnet::xml::tag("FullMFT") << magnet::xml::attr("Src1")
        << getEventSourceName(k.first.first.first, Sim)
        << magnet::xml::attr("Event1") << k.first.first.second
        << magnet::xml::attr("Captures1") << k.first.second
        << magnet::xml::attr("Src2")
        << getEventSourceName(k.second.first.first, Sim)
        << magnet::xml::attr("Event2") << k.second.first.second
        << magnet::xml::attr("Captures2") << k.second.second;
    h.outputHistogram(XML, 1.0 / Sim->units.unitTime());
    XML << magnet::xml::endtag("FullMFT");
  }
}
} // namespace dynamo

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

#pragma once
#include <dynamo/eventtypes.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/math/histogram.hpp>
#include <map>
#include <vector>

namespace dynamo {
class Particle;

using namespace EventTypeTracking;

class OPCollMatrix : public OutputPlugin {
private:
public:
  OPCollMatrix(const dynamo::Simulation *, const magnet::xml::Node &);
  ~OPCollMatrix();

  virtual void initialise();

  virtual void eventUpdate(const Event &, const NEventData &);

  void output(magnet::xml::XmlStream &);

protected:
  void newEvent(const size_t &, const EEventType &, const EventSourceKey &);

  struct InterEventData {
    InterEventData() : count(0), totalTime(0) {}
    unsigned long count;
    double totalTime;
  };

  unsigned long totalCount;

  // EventKet is a pair of EventSourceKey and EEventType
  // It describes the type of event and its source

  //! A key for two events
  typedef std::pair<EventKey, EventKey> InterEventKey;

  //! Counters for properties between two different events
  std::map<InterEventKey, InterEventData> counters;

  // First we track how many times a particle has been captured
  typedef std::pair<size_t, size_t> TotalCaptureStateKey; // Interaction ID and particle ID
  struct CaptureStateData {
    CaptureStateData(double binWidth = 1.0) {}
    double _last_update = 0;
    size_t _state = 0;
    EventKey _last_event = EventKey(EventSourceKey(0, NOSOURCE), NONE);
    double _last_event_time = 0;
  };
  std::map<TotalCaptureStateKey, CaptureStateData>
      _currentCaptureState; // How many captures a particle has

  magnet::math::HistogramWeighted<> _captureStateHistogram;

  // Here we're tracking collision statistics depending on the Event Type/Source
  // and capture state of each particle individually
  typedef std::pair<EventKey, size_t> EventCaptureStateKey;
  struct EventCaptureStateData {
    EventCaptureStateData(double binWidth)
        : MFT(binWidth), rijdotvij(0.01), rijdotdP(0.01), vi2(0.01) {}
    double last_event_time = 0;
    magnet::math::Histogram<> MFT;
    magnet::math::Histogram<> rijdotvij;
    magnet::math::Histogram<> rijdotdP;
    magnet::math::Histogram<> vi2;
    magnet::math::Histogram<> _particle_MFT;
  };
  std::map<EventCaptureStateKey, EventCaptureStateData> _captureCounters;

  typedef std::tuple<EventKey, size_t, size_t> PairEventCaptureStateKey;
  struct PairEventCaptureStateData {
    PairEventCaptureStateData(double binWidth)
        : MFT(binWidth), rijdotvij(0.01), rijdotdP(0.01), vi2(0.01) {}
    double last_event_time = 0;
    magnet::math::Histogram<> MFT;
    magnet::math::Histogram<> rijdotvij;
    magnet::math::Histogram<> rijdotdP;
    magnet::math::Histogram<> vi2;
    magnet::math::Histogram<> _particle_MFT;
  };
  std::map<PairEventCaptureStateKey, PairEventCaptureStateData> _pairCaptureCounters;



  typedef std::pair<EventCaptureStateKey, EventCaptureStateKey> MFTKey;

  std::map<MFTKey, magnet::math::Histogram<>> _fullMFT;

  std::map<EventKey, size_t> initialCounter;

  typedef std::pair<double, EventKey> lastEventData;

  std::vector<lastEventData> lastEvent;
};
} // namespace dynamo

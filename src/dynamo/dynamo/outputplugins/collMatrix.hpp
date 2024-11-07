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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/eventtypes.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <map>
#include <vector>

namespace dynamo {
  class Particle;

  using namespace EventTypeTracking;

  class OPCollMatrix: public OutputPlugin
  {
  private:
  
  public:
    OPCollMatrix(const dynamo::Simulation*, const magnet::xml::Node&);
    ~OPCollMatrix();

    virtual void initialise();

    virtual void eventUpdate(const Event&, const NEventData&);

    void output(magnet::xml::XmlStream &);
  
  protected:
    void newEvent(const size_t&, const EEventType&, const classKey&);
  
    struct counterData
    {
      counterData():count(0), totalTime(0) {}
      unsigned long count;
      double totalTime;
    };
  
    unsigned long totalCount;

    // We create a key for events based on the interaction/system/global/local ID and type (classKey) and EventType
    typedef std::pair<classKey, EEventType> eventKey;

    // And we count time between events for a particular particle, so we need this twice
    typedef std::pair<eventKey, eventKey> counterKey;
  
    std::map<counterKey, counterData> counters;

    //First we track how many times a particle has been captured
    typedef std::pair<size_t, size_t> captureKey; // Interaction ID and particle ID
    std::map<captureKey, size_t> _lastCaptureState; // How many captures a particle has


    //Here we're tracking collision statistics depending on the capture state of
    //a pair of particles
    typedef std::pair<eventKey, size_t> captureEventKey;
    struct captureData {
      captureData():count(0) {}
      unsigned long count;
    };
    std::map<captureEventKey, captureData> _captureCounters;

    std::map<eventKey, size_t> initialCounter;

    typedef std::pair<double, eventKey> lastEventData;

    std::vector<lastEventData> lastEvent; 
  };
}

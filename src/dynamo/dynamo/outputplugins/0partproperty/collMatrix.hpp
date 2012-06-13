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

    virtual void eventUpdate(const IntEvent&, const PairEventData&);

    virtual void eventUpdate(const GlobalEvent&, const NEventData&);

    virtual void eventUpdate(const LocalEvent&, const NEventData&);

    virtual void eventUpdate(const System&, const NEventData&, const double&);

    void output(magnet::xml::XmlStream &);

    //This is fine to replica exchange as the interaction, global and system lookups are done using names
    virtual void changeSystem(OutputPlugin* plug) { std::swap(Sim, static_cast<OPCollMatrix*>(plug)->Sim); }
  
  protected:
    void newEvent(const size_t&, const EEventType&, const classKey&);
  
    struct counterData
    {
      counterData():count(0), initialCount(0), totalTime(0) {}
      unsigned long count;
      size_t initialCount;
      double totalTime;
    };
  
    unsigned long totalCount;

    typedef std::pair<classKey, EEventType> eventKey;

    typedef std::pair<eventKey, eventKey> counterKey;
  
    std::map<counterKey, counterData> counters;
  
    std::map<eventKey, size_t> initialCounter;

    typedef std::pair<double, eventKey> lastEventData;

    std::vector<lastEventData> lastEvent; 
  };
}

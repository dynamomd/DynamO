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

#ifndef OPCollMatrix_H
#define OPCollMatrix_H

#include "../outputplugin.hpp"
#include "../../dynamics/eventtypes.hpp"
#include <map>
#include <vector>
#include "../eventtypetracking.hpp"

using namespace EventTypeTracking;

class CParticle;

class OPCollMatrix: public OutputPlugin
{
private:
  
public:
  OPCollMatrix(const DYNAMO::SimData*, const XMLNode&);
  ~OPCollMatrix();

  virtual void initialise();

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  void output(xmlw::XmlStream &);

  virtual OutputPlugin *Clone() const { return new OPCollMatrix(*this); };

  //This is fine to replica exchange as the interaction, global and system lookups are done using names
  virtual void changeSystem(OutputPlugin* plug) { std::swap(Sim, static_cast<OPCollMatrix*>(plug)->Sim); }
  
 protected:
  void newEvent(const size_t&, const EEventType&, const classKey&);
  
  struct counterData
  {
    counterData():count(0), initialCount(0), totalTime(0) {}
    unsigned long count;
    size_t initialCount;
    Iflt totalTime;
  };
  
  unsigned long totalCount;

  typedef std::pair<classKey, EEventType> eventKey;

  typedef std::pair<eventKey, eventKey> counterKey;
  
  std::map<counterKey, counterData> counters;
  
  std::map<eventKey, size_t> initialCounter;

  typedef std::pair<Iflt, eventKey> lastEventData;

  std::vector<lastEventData> lastEvent; 
};

#endif

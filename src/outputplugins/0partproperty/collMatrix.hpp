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

#ifndef COPCollMatrix_H
#define COPCollMatrix_H

#include "../outputplugin.hpp"
#include "../../dynamics/eventtypes.hpp"
#include <map>
#include <vector>

class CParticle;

class COPCollMatrix: public COutputPlugin
{
 public:
  COPCollMatrix(const DYNAMO::SimData*);
  ~COPCollMatrix();

  virtual void initialise();

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  void output(xmlw::XmlStream &);

  virtual COutputPlugin *Clone() const { return new COPCollMatrix(*this); };

  //This is fine to replica exchange as the interaction, global and system lookups are done using names
  virtual void changeSystem(COutputPlugin* plug) { std::swap(Sim, static_cast<COPCollMatrix*>(plug)->Sim); }
  
  size_t getID(const CInteraction&) const;
  size_t getID(const CGlobal&) const;
  size_t getID(const CSystem&) const;
  std::string getName(const unsigned int&) const;

 protected:
  void newEvent(const unsigned int, const CParticle&, EEventType);

  struct counterData
  {
    counterData():count(0),totalTime(0) {}
    unsigned long count;
    Iflt totalTime;
  };
  
  unsigned long totalCount;
  
  mutable std::map<const std::string, unsigned int> intLookup;
  mutable std::map<const std::string, unsigned int> globLookup;
  mutable std::map<const std::string, unsigned int> sysLookup;

  typedef std::pair<std::pair<unsigned int, EEventType>, std::pair<unsigned int, EEventType> > counterKey;
  typedef std::pair<const counterKey, counterData> counterElement;
  std::map<counterKey, counterData> counters;
  
  typedef std::pair<Iflt, std::pair<unsigned int, EEventType> > p2timeData;
  std::vector<p2timeData> p2time; 
  mutable size_t IDcounter;
};

#endif

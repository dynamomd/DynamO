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

#pragma once
#include "../base/is_base.hpp"
#include "sorters/cbt.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include <vector>

class Particle;
namespace magnet { namespace xml { class Node; } }

class CScheduler: public dynamo::SimBase
{
public:
  CScheduler(dynamo::SimData* const, const char *, CSSorter*);
  
  virtual ~CScheduler() = 0;

  virtual void initialise() = 0;
  
  virtual void fullUpdate(const Particle& part);

  virtual void fullUpdate(const Particle& p1, const Particle& p2);

  void invalidateEvents(const Particle&);

  virtual void addEvents(const Particle&) = 0;

  void sort(const Particle&);

  void popNextEvent();

  void pushEvent(const Particle&, const intPart&);
  
  void stream(const double& dt) {  sorter->stream(dt); }
  
  void runNextEvent();

  virtual void rebuildList() = 0;

  friend xml::XmlStream& operator<<(xml::XmlStream&, const CScheduler&);

  static CScheduler* getClass(const magnet::xml::Node&, dynamo::SimData* const);

  virtual void operator<<(const magnet::xml::Node&) = 0;
  
  void rescaleTimes(const double& scale) { sorter->rescaleTimes(scale); }

  const magnet::ClonePtr<CSSorter>& getSorter() const { return sorter; }

  void rebuildSystemEvents() const;

  void addInteractionEvent(const Particle&, const size_t&) const;

  void addInteractionEventInit(const Particle&, const size_t&) const;

  void addLocalEvent(const Particle&, const size_t&) const;
  
protected:
  mutable magnet::ClonePtr<CSSorter> sorter;
  mutable std::vector<unsigned long> eventCount;
  
  size_t _interactionRejectionCounter;
  size_t _localRejectionCounter;

  virtual void outputXML(xml::XmlStream&) const = 0;
};

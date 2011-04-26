/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "../../base/is_base.hpp"
#include "../ranges/2range.hpp"
#include <magnet/cloneptr.hpp>
#include <string>

class PairEventData;
class IntEvent;
class Species;
class CRange;

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

class Interaction: public DYNAMO::SimBase
{
public:
  Interaction(DYNAMO::SimData*, C2Range*);
  
  virtual ~Interaction() {}

  virtual void initialise(size_t) = 0;

  virtual IntEvent getEvent(const Particle &, 
			     const Particle &) const = 0;

  virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const = 0;

  virtual double maxIntDist() const = 0;  

  //! Returns the internal energy "stored" in this interaction
  virtual double getInternalEnergy() const = 0; 

  virtual double getExcludedVolume(size_t) const = 0;

  virtual Interaction* Clone() const = 0; //{ return new OPBlank(*this); };

  virtual void operator<<(const magnet::xml::Node&) = 0;
  
  friend xml::XmlStream& operator<<(xml::XmlStream&, const Interaction&);
 
  static Interaction* getClass(const magnet::xml::Node&, DYNAMO::SimData*);

  bool isInteraction(const Particle &p1, const Particle &p2) const
  { return range->isInRange(p1,p2); }
  
  bool isInteraction(const IntEvent &) const;

  bool isInteraction(const Species &) const;

  inline void setName(const std::string& tmp) { intName = tmp; }

  inline const std::string& getName() const { return intName; }

  magnet::ClonePtr<C2Range>& getRange();

  const magnet::ClonePtr<C2Range>& getRange() const;

  virtual void checkOverlaps(const Particle&, const Particle&) const = 0;

  virtual double getColourFraction(const Particle&) const { return 0.5; } 

  inline const size_t& getID() const { return ID; }

protected:
  virtual void outputXML(xml::XmlStream& ) const = 0;

  magnet::ClonePtr<C2Range> range;

  std::string intName;
  size_t ID;
};

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

#ifndef CInteraction_H
#define CInteraction_H

#include <string>
#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/2range.hpp"
#include "../../base/is_colormap.hpp"

#define eps2 1e-10 //Min overlap distance

class C2ParticleData;
class CIntEvent;
class XMLNode;
class CSpecies;
class CRange;

namespace xmlw
{
  class XmlStream;
}

class CInteraction: public DYNAMO::SimBase
{
public:
  CInteraction(DYNAMO::SimData*, C2Range*);
  
  virtual ~CInteraction() {}

  virtual void initialise(size_t) = 0;

  virtual CIntEvent getEvent(const CParticle &, 
			     const CParticle &) const = 0;

  virtual void runEvent(const CParticle&, const CParticle&) const = 0;

  virtual Iflt maxIntDist() const = 0;  

  virtual Iflt getInternalEnergy() const = 0; 

  virtual Iflt hardCoreDiam() const = 0;

  virtual void rescaleLengths(Iflt) = 0;

  virtual CInteraction* Clone() const = 0; //{ return new COPBlank(*this); };

  virtual void operator<<(const XMLNode&) = 0;
  
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CInteraction&);
 
  static CInteraction* getClass(const XMLNode&, DYNAMO::SimData*);

  bool isInteraction(const CParticle &p1, const CParticle &p2) const
  { return range->isInRange(p1,p2); }
  
  bool isInteraction(const CIntEvent &) const;

  bool isInteraction(const CSpecies &) const;

  inline void setName(const std::string& tmp) { intName = tmp; }

  inline const std::string& getName() const { return intName; }

  smrtPlugPtr<C2Range>& getRange();

  const smrtPlugPtr<C2Range>& getRange() const;

  virtual void checkOverlaps(const CParticle&, const CParticle&) const = 0;

  virtual Iflt getColourFraction(const CParticle&) const { return 0.5; } 

  inline const size_t& getID() const { return ID; }

  virtual void 
  write_povray_desc(const DYNAMO::RGB&, const size_t&, std::ostream&) const
  {}

  virtual void 
  write_povray_info(std::ostream&) const 
  {}

protected:
  virtual void outputXML(xmlw::XmlStream& ) const = 0;

  smrtPlugPtr<C2Range> range;

  std::string intName;
  size_t ID;
};

#endif

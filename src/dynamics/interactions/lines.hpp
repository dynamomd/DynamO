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

#ifndef CILines_H
#define CILines_H

#include "captures.hpp"

class CILines: public CICapture
{
public:
  CILines(DYNAMO::SimData*, Iflt, Iflt, C2Range*);

  CILines(const XMLNode&, DYNAMO::SimData*);

  void operator<<(const XMLNode&);

  virtual Iflt getInternalEnergy() const { return 0.0; }

  virtual void initialise(size_t);

  virtual Iflt maxIntDist() const;

  virtual Iflt hardCoreDiam() const;

  virtual void rescaleLengths(Iflt);

  virtual CInteraction* Clone() const;
  
  virtual CIntEvent getEvent(const CParticle&, const CParticle&) const;
 
  virtual void runEvent(const CParticle&, const CParticle&, const CIntEvent&) const;
   
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void checkOverlaps(const CParticle&, const CParticle&) const;
 
  virtual bool captureTest(const CParticle&, const CParticle&) const;

  virtual void write_povray_desc(const DYNAMO::RGB&, 
				 const size_t&, std::ostream&) const;

protected:
  Iflt length;
  Iflt l2;
  Iflt e;
};

#endif

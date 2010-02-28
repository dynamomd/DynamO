/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CISoftCore_H
#define CISoftCore_H

#include "captures.hpp"

class CISoftCore: public CISingleCapture
{
public:
  CISoftCore(DYNAMO::SimData*, Iflt, Iflt, C2Range*);

  CISoftCore(const XMLNode&, DYNAMO::SimData*);
  
  void operator<<(const XMLNode&);

  virtual CInteraction* Clone() const;

  virtual Iflt hardCoreDiam() const;

  virtual Iflt maxIntDist() const;

  virtual void rescaleLengths(Iflt);

  virtual void checkOverlaps(const CParticle&, const CParticle&) const;

  virtual bool captureTest(const CParticle&, const CParticle&) const;

  virtual void initialise(size_t);

  virtual CIntEvent getEvent(const CParticle&, const CParticle&) const;
  
  virtual void runEvent(const CParticle&, const CParticle&, const CIntEvent&) const;
  
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual Iflt getInternalEnergy() const { return -(getTotalCaptureCount() * wellDepth); }

protected:
  Iflt diameter,d2;
  Iflt wellDepth;
};

#endif

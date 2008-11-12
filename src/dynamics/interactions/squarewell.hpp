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

#ifndef CISquareWell_H
#define CISquareWell_H

#include "captures.hpp"

class CISquareWell: public CICapture
{
public:
  CISquareWell(const DYNAMO::SimData*, Iflt, Iflt, Iflt, Iflt, C2Range*);

  CISquareWell(const XMLNode&, const DYNAMO::SimData*);
  
  void operator<<(const XMLNode&);

  virtual CInteraction* Clone() const;

  virtual Iflt hardCoreDiam() const;

  virtual Iflt maxIntDist() const;

  virtual void rescaleLengths(Iflt);

  virtual void checkOverlaps(const CParticle&, const CParticle&) const;

  virtual bool captureTest(const CParticle&, const CParticle&) const;

  virtual void initialise(size_t);

  virtual CIntEvent getCollision(const CParticle&, const CParticle&) const;
  
  virtual C2ParticleData runCollision(const CIntEvent&) const;
  
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual Iflt getInternalEnergy() const 
  { return -(getTotalCaptureCount() * wellDepth); }

  virtual void 
  write_povray_desc(const DYNAMO::RGB&, const CRange&, std::ostream&) const;

protected:
  Iflt diameter,d2;
  Iflt lambda, ld2;
  Iflt wellDepth;
  Iflt e;
};

#endif

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

#ifndef CISWSequence_H
#define CISWSequence_H

#include "captures.hpp"
#include <vector>

class ISWSequence: public ISingleCapture
{
public:
  ISWSequence(DYNAMO::SimData*, Iflt, Iflt, Iflt, std::vector<size_t>, C2Range*);

  ISWSequence(const XMLNode&, DYNAMO::SimData*);
  
  void operator<<(const XMLNode&);

  virtual Interaction* Clone() const;

  virtual Iflt hardCoreDiam() const;

  virtual Iflt maxIntDist() const;

  virtual Iflt getInternalEnergy() const;

  virtual void rescaleLengths(Iflt);

  virtual void checkOverlaps(const Particle&, const Particle&) const;

  virtual bool captureTest(const Particle&, const Particle&) const;

  virtual void initialise(size_t);

  virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
  virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const;
  
  virtual void outputXML(xml::XmlStream&) const;

  virtual Iflt getColourFraction(const Particle&) const;

  std::vector<size_t>& getSequence() { return sequence; }

  std::vector<std::vector<Iflt> >& getAlphabet() { return alphabet; }

  virtual void write_povray_desc(const DYNAMO::RGB&, 
				 const size_t&, 
				 std::ostream&) const;

protected:
  Iflt diameter,d2;
  Iflt lambda, ld2;
  Iflt e;

  std::vector<size_t> sequence;
  std::vector<std::vector<Iflt> > alphabet;
};

#endif

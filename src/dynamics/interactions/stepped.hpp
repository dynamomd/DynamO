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

#ifndef CIStepped_H
#define CIStepped_H

#include "captures.hpp"

class IStepped: public IMultiCapture
{
public:
  typedef std::pair<double,double> steppair;

  IStepped(DYNAMO::SimData*, const std::vector<steppair>&,
	    C2Range*);

  IStepped(const XMLNode&, DYNAMO::SimData*);
  
  void operator<<(const XMLNode&);

  virtual Interaction* Clone() const;

  virtual double hardCoreDiam() const;

  virtual double maxIntDist() const;

  virtual void rescaleLengths(double);

  virtual void checkOverlaps(const Particle&, const Particle&) const;

  virtual int captureTest(const Particle&, const Particle&) const;

  virtual void initialise(size_t);

  virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
  virtual void runEvent(const Particle&, const Particle&, 
			const IntEvent&) const;
  
  virtual void outputXML(xml::XmlStream&) const;

  virtual double getInternalEnergy() const;

  virtual void 
  write_povray_desc(const DYNAMO::RGB&, const size_t&, std::ostream&) const;

protected:
  std::vector<steppair> steps;

  std::vector<steppair> runstepdata;
};

#endif

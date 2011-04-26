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

#include "captures.hpp"
#include "../../base/is_simdata.hpp"
#include <vector>

class IStepped: public IMultiCapture
{
public:
  typedef std::pair<double,double> steppair;

  IStepped(DYNAMO::SimData*, const std::vector<steppair>&,
	    C2Range*);

  IStepped(const magnet::xml::Node&, DYNAMO::SimData*);
  
  void operator<<(const magnet::xml::Node&);

  virtual Interaction* Clone() const;

  virtual double getExcludedVolume(size_t) const;

  virtual double maxIntDist() const;

  virtual void checkOverlaps(const Particle&, const Particle&) const;

  virtual int captureTest(const Particle&, const Particle&) const;

  virtual void initialise(size_t);

  virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
  virtual void runEvent(const Particle&, const Particle&, 
			const IntEvent&) const;
  
  virtual void outputXML(xml::XmlStream&) const;

  virtual double getInternalEnergy() const;

protected:
  //!This class is used to track how the length scale changes in the system
  magnet::thread::RefPtr<Property> _unitLength;
  //!This class is used to track how the energy scale changes in the system
  magnet::thread::RefPtr<Property> _unitEnergy;

  std::vector<steppair> steps;
};

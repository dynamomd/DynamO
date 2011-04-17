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

#include "interaction.hpp"
#include "../../base/is_simdata.hpp"

class IHardSphere: public Interaction
{
public:
  template<class T1, class T2>
  IHardSphere(DYNAMO::SimData* tmp, T1 diameter, T2 e, C2Range* nR):
  Interaction(tmp, nR),
  _diameter(Sim->_properties.getProperty
	    (diameter, Property::Units::Length())),
  _e(Sim->_properties.getProperty
     (e, Property::Units::Dimensionless()))
  {}

  IHardSphere(const magnet::xml::Node&, DYNAMO::SimData*);

  void operator<<(const magnet::xml::Node&);

  virtual double getInternalEnergy() const { return 0.0; }

  virtual void initialise(size_t);

  virtual double maxIntDist() const;

  virtual double hardCoreDiam() const;

  virtual void rescaleLengths(double) {}

  virtual Interaction* Clone() const;
  
  virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
  virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const;
   
  virtual void outputXML(xml::XmlStream&) const;

  virtual void checkOverlaps(const Particle&, const Particle&) const;

protected:
  magnet::thread::RefPtr<Property> _diameter;
  magnet::thread::RefPtr<Property> _e;
};

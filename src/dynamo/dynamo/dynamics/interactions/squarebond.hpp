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

#include "interaction.hpp"
#include "../../base/is_simdata.hpp"

class ISquareBond: public Interaction
{
public:
  template<class T1, class T2, class T3>
  ISquareBond(dynamo::SimData* tmp, T1 d, T2 l, T3 e, C2Range* nR):
    Interaction(tmp, nR),
    _diameter(Sim->_properties.getProperty
	      (d, Property::Units::Length())),
    _lambda(Sim->_properties.getProperty
	    (l, Property::Units::Dimensionless())),
    _e(Sim->_properties.getProperty
       (e, Property::Units::Dimensionless()))
  {}

  ISquareBond(const magnet::xml::Node&, dynamo::SimData*);

  void operator<<(const magnet::xml::Node&);

  virtual Interaction* Clone() const;

  virtual double getExcludedVolume(size_t) const 
  { M_throw() << "Bonds don't have excluded volumes! They shouldn't be used as the defining interaction for a species."; }

  virtual double maxIntDist() const;

  virtual double getCaptureEnergy() const;

  virtual void initialise(size_t);

  virtual bool captureTest(const Particle&, const Particle&) const;

  virtual void checkOverlaps(const Particle&, const Particle&) const;

  virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
  virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const;
    
  virtual void outputXML(magnet::xml::XmlStream&) const;

  virtual double getInternalEnergy() const { return 0.0; }

protected:

  magnet::thread::RefPtr<Property> _diameter;
  magnet::thread::RefPtr<Property> _lambda;
  magnet::thread::RefPtr<Property> _e;
};

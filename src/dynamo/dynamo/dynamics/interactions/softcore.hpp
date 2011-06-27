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

#include "captures.hpp"
#include "representations/spherical.hpp"
#include "../../base/is_simdata.hpp"

class ISoftCore: public ISingleCapture, public Interaction, public SphericalRepresentation
{
public:
  template<class T1, class T2>
  ISoftCore(dynamo::SimData* tmp, T1 d, T2 wd, C2Range* nR):
    Interaction(tmp,nR),
    _diameter(Sim->_properties.getProperty
	      (d, Property::Units::Length())),
    _wellDepth(Sim->_properties.getProperty
	       (wd, Property::Units::Energy()))
  {}
  
  ISoftCore(const magnet::xml::Node&, dynamo::SimData*);
  
  virtual size_t spheresPerParticle() const { return 1; }
  virtual double getDiameter(size_t ID, size_t subID) const;
  virtual Vector getPosition(size_t ID, size_t subID) const;

  void operator<<(const magnet::xml::Node&);

  virtual Interaction* Clone() const;

  virtual double getExcludedVolume(size_t) const { return 0; }

  virtual double maxIntDist() const;

  virtual void checkOverlaps(const Particle&, const Particle&) const;

  virtual bool captureTest(const Particle&, const Particle&) const;

  virtual void initialise(size_t);

  virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
  virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const;
  
  virtual void outputXML(magnet::xml::XmlStream&) const;

  virtual double getInternalEnergy() const;

protected:
  magnet::thread::RefPtr<Property> _diameter;
  magnet::thread::RefPtr<Property> _wellDepth;
};

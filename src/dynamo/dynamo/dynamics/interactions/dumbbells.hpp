/*  dynamo:- Event driven molecular dynamics simulator 
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

#pragma once

#include <dynamo/dynamics/interactions/captures.hpp>
#include <dynamo/dynamics/interactions/interaction.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/interactions/glyphrepresentation.hpp>

namespace dynamo {
  class IDumbbells: public ISingleCapture, public Interaction, public GlyphRepresentation
  {
  public:
    template<class T1, class T2, class T3>
    IDumbbells(dynamo::SimData* tmp, T1 l, T2 e, T3 d, C2Range* nR, 
	       std::string name):
      Interaction(tmp, nR),
      _length(Sim->_properties.getProperty
	      (l, Property::Units::Length())),
      _diameter(Sim->_properties.getProperty
		(d, Property::Units::Length())),
      _e(Sim->_properties.getProperty
	 (e, Property::Units::Dimensionless()))
    {
      intName = name;
    }

    virtual size_t glyphsPerParticle() const { return 2; }
    virtual double getGlyphDiameter(size_t ID, size_t subID) const;
    virtual Vector getGlyphPosition(size_t ID, size_t subID) const;

    IDumbbells(const magnet::xml::Node&, dynamo::SimData*);

    void operator<<(const magnet::xml::Node&);
  
    double getDiameter() const { return _diameter->getMaxValue(); }

    double getLength() const { return _length->getMaxValue(); }
  
    virtual double getInternalEnergy() const { return 0; }

    virtual void initialise(size_t);

    virtual double maxIntDist() const;

    virtual double getExcludedVolume(size_t) const;

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
    virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const;
   
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual void checkOverlaps(const Particle&, const Particle&) const;
 
    virtual bool captureTest(const Particle&, const Particle&) const;

  protected:
    shared_ptr<Property> _length;
    shared_ptr<Property> _diameter;
    shared_ptr<Property> _e;
  };
}

/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/interactions/captures.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/interactions/glyphrepresentation.hpp>

namespace dynamo {
  class ILines: public ISingleCapture, public GlyphRepresentation
  {
  public:
    template<class T1, class T2>
    ILines(dynamo::Simulation* tmp, T1 l, T2 e, C2Range* nR, std::string name):
      ISingleCapture(tmp, nR),
      _length(Sim->_properties.getProperty
	      (l, Property::Units::Length())),
      _e(Sim->_properties.getProperty
	 (e, Property::Units::Dimensionless()))
    {
      intName = name;
    }

    virtual size_t glyphsPerParticle() const { return 1; }
    virtual Vector getGlyphSize(size_t ID, size_t subID) const;
    virtual Vector getGlyphPosition(size_t ID, size_t subID) const;
    virtual GLYPH_TYPE getDefaultGlyphType() const { return LINE_GLYPH; }

    ILines(const magnet::xml::Node&, dynamo::Simulation*);

    void operator<<(const magnet::xml::Node&);

    virtual double getInternalEnergy() const { return 0; }

    virtual void initialise(size_t);

    virtual double maxIntDist() const;

    virtual double getExcludedVolume(size_t) const { return 0; }

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
    virtual void runEvent(Particle&, Particle&, const IntEvent&) const;
   
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;
 
    virtual bool captureTest(const Particle&, const Particle&) const;

  protected:
    shared_ptr<Property> _length;
    shared_ptr<Property> _e;
  };
}

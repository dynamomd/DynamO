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
#include <dynamo/interactions/glyphrepresentation.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
  class ISoftCore: public ISingleCapture, public GlyphRepresentation
  {
  public:
    template<class T1, class T2>
    ISoftCore(dynamo::Simulation* tmp, T1 d, T2 wd, IDPairRange* nR, std::string name):
      ISingleCapture(tmp,nR),
      _diameter(Sim->_properties.getProperty
		(d, Property::Units::Length())),
      _wellDepth(Sim->_properties.getProperty
		 (wd, Property::Units::Energy()))
    { intName = name; }
  
    ISoftCore(const magnet::xml::Node&, dynamo::Simulation*);
  
    virtual size_t glyphsPerParticle() const { return 1; }
    virtual Vector getGlyphSize(size_t ID, size_t subID) const;
    virtual Vector getGlyphPosition(size_t ID, size_t subID) const;

    void operator<<(const magnet::xml::Node&);

    virtual double getExcludedVolume(size_t) const { return 0; }

    virtual double maxIntDist() const;

    virtual bool captureTest(const Particle&, const Particle&) const;

    virtual void initialise(size_t);

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
    virtual void runEvent(Particle&, Particle&, const IntEvent&) const;
  
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual double getInternalEnergy() const;

    virtual double getInternalEnergy(const Particle&, const Particle&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

    virtual size_t validateState(bool textoutput = true, size_t max_reports = std::numeric_limits<size_t>::max()) const;

  protected:
    shared_ptr<Property> _diameter;
    shared_ptr<Property> _wellDepth;
  };
}


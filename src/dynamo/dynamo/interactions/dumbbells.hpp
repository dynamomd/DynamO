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

namespace dynamo {
  class IDumbbells: public ICapture
  {
  public:
    template<class T1, class T2, class T3, class T4, class T5>
    IDumbbells(dynamo::Simulation* tmp, T1 LA, T2 LB, T3 diamA, T4 diamB, T5 e, IDPairRange* nR, std::string name):
      ICapture(tmp, nR),
      _diamA(Sim->_properties.getProperty(diamA, Property::Units::Length())),
      _diamB(Sim->_properties.getProperty(diamB, Property::Units::Length())),
      _LA(Sim->_properties.getProperty(LA, Property::Units::Length())),
      _LB(Sim->_properties.getProperty(LB, Property::Units::Length())),
      _e(Sim->_properties.getProperty(e, Property::Units::Dimensionless()))
    {
      intName = name;
    }

    virtual std::array<double, 4> getGlyphSize(size_t ID) const;

    virtual GLYPH_TYPE getDefaultGlyphType() const { return DUMBBELL_GLYPH; }

    virtual double getExcludedVolume(size_t ID) const;

    IDumbbells(const magnet::xml::Node&, dynamo::Simulation*);

    void operator<<(const magnet::xml::Node&);

    virtual void initialise(size_t);

    virtual double maxIntDist() const;

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
    virtual PairEventData runEvent(Particle&, Particle&, const IntEvent&);
   
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual size_t captureTest(const Particle&, const Particle&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

  protected:
    shared_ptr<Property> _diamA;
    shared_ptr<Property> _diamB;
    shared_ptr<Property> _LA;
    shared_ptr<Property> _LB;
    shared_ptr<Property> _e;
  };
}

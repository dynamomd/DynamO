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

#include <dynamo/interactions/interaction.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/interactions/glyphrepresentation.hpp>

namespace dynamo {
  class IHardSphere: public GlyphRepresentation, public Interaction
  {
  public:
    template<class T1, class T2>
    IHardSphere(dynamo::Simulation* tmp, T1 d, T2 e, IDPairRange* nR, 
		std::string name):
      Interaction(tmp, nR),
      _diameter(Sim->_properties.getProperty
		(d, Property::Units::Length())),
      _e(Sim->_properties.getProperty
	 (e, Property::Units::Dimensionless()))
    {
      intName = name;
    }

    IHardSphere(const magnet::xml::Node&, dynamo::Simulation*);

    void operator<<(const magnet::xml::Node&);

    virtual Vector getGlyphSize(size_t ID) const;

    virtual void initialise(size_t);

    virtual double maxIntDist() const;

    virtual double getExcludedVolume(size_t) const;

    virtual void rescaleLengths(double) {}

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
    virtual void runEvent(Particle&, Particle&, const IntEvent&);
   
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

    void outputData(magnet::xml::XmlStream& XML) const;

  protected:
    shared_ptr<Property> _diameter;
    shared_ptr<Property> _e;
    
    mutable size_t _complete_events;
    mutable size_t _post_event_overlap;
    mutable double _accum_overlap_magnitude;
    mutable size_t _overlapped_tests;
  };
}

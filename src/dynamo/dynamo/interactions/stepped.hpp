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
#include <dynamo/interactions/potentials/potential.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/eventtypes.hpp>
#include <vector>

namespace dynamo {
  class IStepped: public ICapture, public GlyphRepresentation
  {
  public:
    template<class T1, class T2>
    IStepped(dynamo::Simulation* tmp, shared_ptr<Potential> potential, IDPairRange* nR, std::string name, T1 length, T2 energy):
      ICapture(tmp,nR),
      _lengthScale(Sim->_properties.getProperty(length, Property::Units::Length())),
      _energyScale(Sim->_properties.getProperty(energy, Property::Units::Energy())),
      _potential(potential)
    {
      intName = name;
    }

    IStepped(const magnet::xml::Node&, dynamo::Simulation*);
  
    void operator<<(const magnet::xml::Node&);

    virtual Vector getGlyphSize(size_t ID) const;

    virtual double getExcludedVolume(size_t) const;

    virtual double maxIntDist() const;

    virtual size_t captureTest(const Particle&, const Particle&) const;

    virtual void initialise(size_t);

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
    virtual void runEvent(Particle&, Particle&, const IntEvent&);
  
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual double getInternalEnergy() const;

    virtual double getInternalEnergy(const Particle&, const Particle&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

    virtual void outputData(magnet::xml::XmlStream&) const;

  protected:
    //!This class is used to track how the length scale changes in the system
    shared_ptr<Property> _lengthScale;
    //!This class is used to track how the energy scale changes in the system
    shared_ptr<Property> _energyScale;

    shared_ptr<Potential> _potential;
    
    struct EdgeData {
      EdgeData(): counter(0), rdotv_sum(0) {}
      size_t counter;
      double rdotv_sum;
    };
    std::map<std::pair<size_t, EEventType>, EdgeData> _edgedata;
  };
}

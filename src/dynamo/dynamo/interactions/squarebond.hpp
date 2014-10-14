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

namespace dynamo {
  class ISquareBond: public Interaction
  {
  public:
    template<class T1, class T2, class T3>
    ISquareBond(dynamo::Simulation* tmp, T1 d, T2 l, T3 e, IDPairRange* nR,
		std::string name):
      Interaction(tmp, nR),
      _diameter(Sim->_properties.getProperty
		(d, Property::Units::Length())),
      _lambda(Sim->_properties.getProperty
	      (l, Property::Units::Dimensionless())),
      _e(Sim->_properties.getProperty
	 (e, Property::Units::Dimensionless()))
    {
      intName = name;
    }

    ISquareBond(const magnet::xml::Node&, dynamo::Simulation*);

    void operator<<(const magnet::xml::Node&);

    virtual double getExcludedVolume(size_t) const 
    { M_throw() << "Bonds don't have excluded volumes! They shouldn't be used as the defining interaction for a species."; }

    virtual double maxIntDist() const;

    virtual double getCaptureEnergy() const;

    virtual bool captureTest(const Particle&, const Particle&) const;

    virtual Event getEvent(const Particle&, const Particle&) const;
  
    virtual PairEventData runEvent(Particle&, Particle&, Event);
    
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

    virtual size_t validateState(bool textoutput = true, size_t max_reports = std::numeric_limits<size_t>::max()) const;
 
  protected:

    shared_ptr<Property> _diameter;
    shared_ptr<Property> _lambda;
    shared_ptr<Property> _e;
  };
}

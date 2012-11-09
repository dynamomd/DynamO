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

namespace dynamo {
  class INull: public Interaction
  {
  public:
    INull(dynamo::Simulation*, C2Range*, std::string);

    INull(const magnet::xml::Node&, dynamo::Simulation*);

    void operator<<(const magnet::xml::Node&);

    virtual double getInternalEnergy() const { return 0; }

    virtual void initialise(size_t);

    virtual double maxIntDist() const { return 0; }

    virtual double getExcludedVolume(size_t) const { return 0; }

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
 
    virtual void runEvent(Particle&, Particle&, const IntEvent&) const;
   
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const { return false; }
  };
}


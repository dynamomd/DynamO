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
#include <dynamo/schedulers/scheduler.hpp>

namespace dynamo {
  class SCEntry;

  class SComplex: public Scheduler
  {
  public:
    SComplex(const magnet::xml::Node&, dynamo::Simulation* const);

    SComplex(dynamo::Simulation* const, FEL*);

    virtual void initialise();

    virtual void operator<<(const magnet::xml::Node&);

    virtual std::auto_ptr<IDRange> getParticleNeighbours(const Particle&) const  { M_throw() << "Unimplemented"; }
    virtual std::auto_ptr<IDRange> getParticleNeighbours(const Vector&) const { M_throw() << "Unimplemented"; }
    virtual std::auto_ptr<IDRange> getParticleLocals(const Particle&) const { M_throw() << "Unimplemented"; }
    
  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    std::vector<shared_ptr<SCEntry> > entries;
  };
}

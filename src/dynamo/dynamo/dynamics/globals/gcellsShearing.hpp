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
#include <dynamo/dynamics/globals/gcellsmorton.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>

namespace dynamo {
  class GCellsShearing: public GCells
  {
  public:
    GCellsShearing(const magnet::xml::Node&, dynamo::SimData*);
  
    GCellsShearing(dynamo::SimData*, const std::string&);
  
    virtual ~GCellsShearing() {}

    virtual void initialise(size_t);
  
    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(const Particle&, const double) const;

    virtual void getParticleNeighbourhood(const Particle&,
					  const nbHoodFunc&) const;

    virtual void getParticleNeighbourhood(const Vector&, 
					  const nbHoodFunc2&) const;

    void getExtraLEParticleNeighbourhood(const Particle& part,
					 const nbHoodFunc& func) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}

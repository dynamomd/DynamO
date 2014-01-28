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
#include <dynamo/globals/cells.hpp>
#include <dynamo/ranges/IDRange.hpp>

namespace dynamo {
  class GCellsShearing: public GCells
  {
  public:
    GCellsShearing(const magnet::xml::Node&, dynamo::Simulation*);
  
    GCellsShearing(dynamo::Simulation*, const std::string&);
  
    virtual ~GCellsShearing() {}

    virtual void initialise(size_t);
  
    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(Particle&, const double) const;

  protected:
    void getParticleNeighbours(const std::array<size_t, 3>&, std::vector<size_t>&) const;
    void getAdditionalLEParticleNeighbourhood(const Particle&, std::vector<size_t>&) const;
    void getAdditionalLEParticleNeighbourhood(std::array<size_t, 3>, std::vector<size_t>&) const;
  };
}

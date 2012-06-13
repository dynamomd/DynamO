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

#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <boost/foreach.hpp>
#include <vector>

namespace magnet {    
  namespace xml { 
    class Node; 
  }
}

namespace dynamo {
  class Particle;
  class BoundaryCondition;
  class Particle;
  class NEventData;
  class PairEventData;
  class ParticleEventData;

  class Dynamics: public dynamo::SimBase
  {
  public:
    //Constructors
    Dynamics(dynamo::SimData*);
    
  protected:
    Dynamics(const Dynamics &dyn);
  };
}

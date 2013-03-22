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
#include <dynamo/interactions/potentials/potential.hpp>
#include <limits>

namespace dynamo {
    /*! \brief The Lennard-Jones potential.
    
      This class implements an automatically stepped potential where
      each step is specified through one of many algorithms.
   */
  class PotentialLennardJones :public Potential {
  public:
    PotentialLennardJones(const magnet::xml::Node& XML) {
      operator<<(XML);
    }
      
    virtual std::size_t steps() const {
      return std::numeric_limits<size_t>::max();
    }

    virtual void operator<<(const magnet::xml::Node&);
  
    virtual double hard_core_diameter() const {
      return 0;
    }

    virtual double render_diameter() const {
      return 0;
    }
    
  protected:
    virtual void calculateNextStep() const;

    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}

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
#include <dynamo/BC/BC.hpp>

namespace dynamo {
  /*! \brief An infinite system boundary condition
   * 
   * Performs no rounding at the simulation boundaries. This is useful
   * for isolated polymer simulations but you must remember that
   * positions can overflow eventually.
   */
  class BCNone: virtual public BoundaryCondition
  {
  public:
    BCNone(const dynamo::Simulation*);
  
    virtual ~BCNone();
    
    virtual void applyBC(Vector&)const;

    virtual void applyBC(Vector&, Vector&) const;

    virtual void applyBC(Vector&, const double&) const;

    virtual void update(const double&);

    virtual void outputXML(magnet::xml::XmlStream &XML) const;

    virtual void operator<<(const magnet::xml::Node&);

    virtual void rounding(Vector &) const;
  };
}

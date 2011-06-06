/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
#include "BC.hpp"

class CLEBC 
{ 
public: 
  static inline double shearRate() { return 1; }
};

/*! \brief A simple rectangular Lees-Edwards simple shear boundary
 * condition.
 *
 * See the square version for more details (BCSquareLeesEdwards)
 */
class BCLeesEdwards: virtual public BoundaryCondition, public CLEBC
{
 public:
  BCLeesEdwards(const dynamo::SimData*);

  BCLeesEdwards(const magnet::xml::Node&, const dynamo::SimData*);

  virtual void outputXML(xml::XmlStream&) const;

  virtual void operator<<(const magnet::xml::Node&);

  virtual BoundaryCondition* Clone () const;

  virtual void applyBC(Vector&) const; 

  virtual void applyBC(Vector&, Vector&) const;

  virtual void applyBC(Vector&, const double& dt) const;

  virtual void update(const double&);

 protected:  
  /*! \brief The current offset of the simulation boundary conditions.*/
  double dxd;
};

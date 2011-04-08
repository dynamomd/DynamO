/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include "../../base/is_base.hpp"

class CLEBC {};

/*! \brief A simple rectangular Lees-Edwards simple shear boundary
 * condition.
 *
 * See the square version for more details (BCSquareLeesEdwards)
 */
class BCRectangularLeesEdwards: virtual public BoundaryCondition, public CLEBC
{
 public:
  BCRectangularLeesEdwards(const DYNAMO::SimData*);

  BCRectangularLeesEdwards(const magnet::xml::Node&, const DYNAMO::SimData*);

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


/*! \brief A simple square Lees-Edwards simple shear boundary condition.
 * 
 * This class implements the sliding brick boundary condtion. In this
 * the simulation images above and below the primary image are set in
 * motion. This affects the particle velocities and positions on a
 * transition of the boundary.
 *
 * See BoundaryCondition for a general description of the member functions.
 */
class BCSquareLeesEdwards: virtual public BoundaryCondition, public CLEBC
{
 public:
  BCSquareLeesEdwards(const DYNAMO::SimData*);

  BCSquareLeesEdwards(const magnet::xml::Node&, const DYNAMO::SimData*);

  virtual void outputXML(xml::XmlStream&) const;

  virtual void operator<<(const magnet::xml::Node&);

  virtual BoundaryCondition* Clone () const;

  virtual void applyBC(Vector&) const; 

  virtual void applyBC(Vector&, Vector&) const;

  virtual void applyBC(Vector&, const double&) const;

  virtual void update(const double&);

 protected:  
  /*! \brief The current offset of the simulation boundary conditions.*/
  double dxd;
};

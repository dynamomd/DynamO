/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef BCNull_H
#define BCNull_H

#include "BC.hpp"
#include "../../base/is_base.hpp"

/*! \brief An infinite system boundary condition
 * 
 * Performs no rounding at the simulation boundaries. This is useful
 * for isolated polymer simulations but you must remember that
 * positions can overflow eventually.
 */
class CNullBC: virtual public CBC
{
 public:
  CNullBC(const DYNAMO::SimData*);
  
  virtual ~CNullBC();
    
  virtual void setPBC(CVector<> &)const;

  virtual void setPBC(CVector<> &, CVector<> &) const;

  virtual void setPBC(CVector<> &, Iflt) const;

  virtual void update(const Iflt&);

  virtual void outputXML(xmlw::XmlStream &XML) const;

  virtual void operator<<(const XMLNode&);

  virtual CBC* Clone () const;

  virtual void rounding(CVector<>&) const;
};

#endif

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

#ifndef PBC_H
#define PBC_H

#include "../../base/is_base.hpp"
#include "BC.hpp"

/*! \brief A simple cubic/square periodic boundary condition.
 * 
 * See the CBC base class for member descriptions.
 */
class CSPBC: public CBC
{
public:
  CSPBC(const DYNAMO::SimData*);

  virtual void applyBC(Vector &) const;

  virtual void applyBC(Vector &, Vector &) const;

  virtual void applyBC(Vector &, const Iflt&) const;

  inline virtual void outputXML(xmlw::XmlStream &) const;  
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

/*! \brief A simple rectangular periodic boundary condition.
 * 
 * See the CBC base class for member descriptions.
 */
class CRPBC: public CBC
{
public:
  CRPBC(const DYNAMO::SimData*);

  virtual void applyBC(Vector &) const;
  
  virtual void applyBC(Vector &, Vector &) const;

  virtual void applyBC(Vector &, const Iflt&) const;

  virtual void outputXML(xmlw::XmlStream&) const;
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

/*! \brief This class ignores the x direction but is periodic in others.
 *
 * Used to check that a system bounded by walls in the x direction has
 * no leaks as these are not rounded and would show up in animations
 * or inspections.
 */
class CRNoXPBC: public CBC
{
public:
  CRNoXPBC(const DYNAMO::SimData*);

  virtual void applyBC(Vector & pos) const;
  
  virtual void applyBC(Vector & pos, Vector &) const;
  
  virtual void applyBC(Vector  &pos, const Iflt&) const;

  virtual void outputXML(xmlw::XmlStream&) const;
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

#endif

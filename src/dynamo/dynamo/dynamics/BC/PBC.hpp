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

/*! \brief A simple rectangular periodic boundary condition.
 * 
 * See the BoundaryCondition base class for member descriptions.
 */
class BCPeriodic: public BoundaryCondition
{
public:
  BCPeriodic(const dynamo::SimData*);

  virtual void applyBC(Vector &) const;
  
  virtual void applyBC(Vector &, Vector &) const;

  virtual void applyBC(Vector &, const double&) const;

  virtual void outputXML(xml::XmlStream&) const;
  virtual void operator<<(const magnet::xml::Node&);
  virtual BoundaryCondition* Clone () const;
};

/*! \brief This class ignores the x direction but is periodic in others.
 *
 * Used to check that a system bounded by walls in the x direction has
 * no leaks as these are not rounded and would show up in animations
 * or inspections.
 */
class BCPeriodicExceptX: public BoundaryCondition
{
public:
  BCPeriodicExceptX(const dynamo::SimData*);

  virtual void applyBC(Vector& pos) const;
  
  virtual void applyBC(Vector& pos, Vector&) const;
  
  virtual void applyBC(Vector& pos, const double&) const;

  virtual void outputXML(xml::XmlStream&) const;
  virtual void operator<<(const magnet::xml::Node&);
  virtual BoundaryCondition* Clone () const;
};

/*! \brief This class ignores all directions but is periodic in
    the x.
 *
 * Used to check that a system bounded by walls in the x direction has
 * no leaks as these are not rounded and would show up in animations
 * or inspections.
 */
class BCPeriodicXOnly: public BoundaryCondition
{
public:
  BCPeriodicXOnly(const dynamo::SimData*);

  virtual void applyBC(Vector& pos) const;
  
  virtual void applyBC(Vector& pos, Vector&) const;
  
  virtual void applyBC(Vector& pos, const double&) const;

  virtual void outputXML(xml::XmlStream&) const;
  virtual void operator<<(const magnet::xml::Node&);
  virtual BoundaryCondition* Clone () const;
};

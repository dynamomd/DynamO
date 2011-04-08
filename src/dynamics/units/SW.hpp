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
#include "units.hpp"

/*! \brief For running a simulation with a distinct energy scale.
 *
 * Mathematically there is no reason for this class and a simulation
 * can be performed using hard sphere units; however, this class is
 * useful for debugging with a certain energy scale. The unit of this
 * scale can then be set to one by adjusting the time scale, which is
 * the function of this class.
 *
 * Although this class is called USquareWell, this class essentially
 * supports any system where there is an inherent energy scale you
 * want the simulation to run in. This usually means the unit energy
 * is set equal to 1.
 */
class USquareWell: public Units
{
 public:
  USquareWell(const DYNAMO::SimData*); 

  USquareWell(double, double, const DYNAMO::SimData*);

  USquareWell(const magnet::xml::Node&, const DYNAMO::SimData*);

  virtual ~USquareWell();

  virtual double unitLength() const;

  virtual void setUnitLength(double);

  virtual double unitTime() const;
  
  virtual Units* Clone() const;
  
  virtual void operator<<(const magnet::xml::Node&);

  virtual void rescaleLength(double);

 protected:
  virtual void outputXML(xml::XmlStream &) const;

  double UnitOfEnergy;
  double UnitOfLength;
};

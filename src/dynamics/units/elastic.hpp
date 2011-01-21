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

#ifndef Elas_Units_H
#define Elas_Units_H

#include "units.hpp"

/*! \brief For running a simulation in hard sphere units.
 *
 * Hard sphere units take a length and mass scale from the diameter of
 * one of the species of particles. The unit of time is arbitrary as
 * the hard sphere system scales trivially with the temperature, so it
 * is typically set such that the temperature is one (This is not done
 * by this class, it will happily work at any temperature. You can
 * scale the temperature to 1 using dynamod).
 */
class UHardSphere: public Units
{
 public:
  UHardSphere(const DYNAMO::SimData*);
  
  UHardSphere(double, const DYNAMO::SimData*);

  UHardSphere(const XMLNode&, const DYNAMO::SimData*);

  virtual ~UHardSphere();

  virtual double unitLength() const;

  virtual void setUnitLength(double);

  virtual double unitTime() const;

  virtual void rescaleLength(double);
  
  virtual Units* Clone() const;
  
  virtual void operator<<(const XMLNode&);

 protected:
  virtual void outputXML(xml::XmlStream&) const;

  double UnitOfLength;
};

#endif

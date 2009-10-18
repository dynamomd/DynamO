/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef SW_Units_H
#define SW_Units_H

#include "units.hpp"

/*! \brief For running a simulation with a distinct energy scale.
 *
 * Scientifically there is no reason for this class and a simulation
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

  USquareWell(Iflt, Iflt, const DYNAMO::SimData*);

  USquareWell(const XMLNode&, const DYNAMO::SimData*);

  virtual ~USquareWell();

  virtual Iflt unitLength() const;

  virtual void setUnitLength(Iflt);

  virtual Iflt unitTime() const;
  
  virtual Units* Clone() const;
  
  virtual void operator<<(const XMLNode&);

  virtual void rescaleLength(Iflt);

 protected:
  virtual void outputXML(xmlw::XmlStream &) const;

  Iflt UnitOfEnergy;
  Iflt UnitOfLength;
};

#endif

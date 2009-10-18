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

#ifndef Units_H
#define Units_H

#include "../../base/is_base.hpp"
#include "../../base/constants.hpp"

class XMLNode;
namespace DYNAMO
{
  class SimData;
}
namespace xmlw
{
  class  XmlStream;
}

/*! \brief The class used to convert in and out of simulation units.
 *
 * The derived classes of Units control the units that the simulation
 * is performed in. All classes scale the system such that its largest
 * box length is unity to allow an optimisation in the square
 * BoundaryCondition.
 *
 * Another reason to set the units is to counteract the above scaling
 * for debugging reasons. Changing the length means you can either
 * preserve the unit time or the unit energy but not both. UShear and
 * USquareWell preserve the time and energy respectively.
 *
 * This class has a particular design specification in that it is
 * initialised on construction so that other classes may directly
 * start conversions on the loading of configurations.
 */
class Units: public DYNAMO::SimBase_const
{
 public:
  Units(const DYNAMO::SimData* const tmp): 
    SimBase_const(tmp, "Units",IC_blue)
    {};

  virtual ~Units() {};

  virtual Units* Clone() const = 0;

  virtual Iflt unitTime() const = 0;
  virtual Iflt unitLength() const = 0;
  virtual void setUnitLength(Iflt) = 0;

  virtual Iflt unitMass() const
    { return 1.0; }
  virtual Iflt unitk() const
    { return 1.0; }

  inline Iflt unitVelocity() const
    { return unitLength() / unitTime(); }
  inline Iflt unitEnergy() const
    { return unitMass() * unitVelocity() * unitVelocity(); }
  inline Iflt unitArea() const
    { return unitLength() * unitLength(); }
  inline Iflt unitVolume() const
    { return unitLength() * unitLength() * unitLength(); }
  inline Iflt unitMomentum() const
    { return unitMass() * unitVelocity();}

  //Some dimensions of some properties
  inline Iflt unitDiffusion() const
  { return unitArea() / unitTime(); }
  inline Iflt unitMutualDiffusion() const
  { return unitMass() * unitTime()/ unitVolume(); }
  inline Iflt unitThermalCond() const
  { return unitk()/(unitLength() * unitTime()); }
  inline Iflt unitThermalDiffusion() const
  { return unitMass()/(unitLength() * unitTime()); }  
  inline Iflt unitViscosity() const
  { return 1.0 /(unitLength() * unitTime()); }
  inline Iflt unitPressure() const
  { return unitMass() / (unitLength()*unitTime()*unitTime()); }

  virtual void rescaleLength(Iflt) = 0;

  Iflt simVolume() const;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const Units&);
  
  virtual void operator<<(const XMLNode &) = 0;
  
  static Units* loadUnits(const XMLNode&, const DYNAMO::SimData*);
  
 protected:
  virtual void outputXML(xmlw::XmlStream &) const = 0;
};

#endif

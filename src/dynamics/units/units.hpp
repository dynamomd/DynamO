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

#ifndef Units_H
#define Units_H

#include "../../base/is_base.hpp"
#include "../../base/constants.hpp"

struct XMLNode;
namespace DYNAMO
{
  class SimData;
}
namespace xml
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

  /*! \brief Pure virtual function returning the simulation unit of time */
  virtual double unitTime() const = 0;

  /*! \brief Pure virtual function returning the simulation unit of length */
  virtual double unitLength() const = 0;

  /*! \brief Pure virtual function returning the method to change the unit length */
  virtual void setUnitLength(double) = 0;

  /*! \brief Overridable simulation unit of mass */
  virtual double unitMass() const
    { return 1.0; }
  
  /*! \brief Overridable Boltzmanns constant */
  virtual double unitk() const
    { return 1.0; }

  /*! Helper function to generate the unit of velocity*/
  inline double unitVelocity() const
    { return unitLength() / unitTime(); }
  inline double unitAcceleration() const
  { return unitLength() / (unitTime() * unitTime()); }
  /*! Helper function to generate the unit of energy*/
  inline double unitEnergy() const
    { return unitMass() * unitVelocity() * unitVelocity(); }
  /*! Helper function to generate the unit of area*/
  inline double unitArea() const
    { return unitLength() * unitLength(); }
  /*! Helper function to generate the unit of volume*/
  inline double unitVolume() const
    { return unitLength() * unitLength() * unitLength(); }
  /*! Helper function to generate the unit of momentum*/
  inline double unitMomentum() const
    { return unitMass() * unitVelocity();}

  //Some dimensions of some properties
  /*! Helper function to generate the units of diffusion outputted by
      the output plugins (see OPMSD/OPMSDCorrelator)*/
  inline double unitDiffusion() const
  { return unitArea() / unitTime(); }

  /*! Helper function to generate the units of mutal diffusion
      outputted by the output plugins (see OPMutualDiffusionE/OPMutualDiffusionGK)*/
  inline double unitMutualDiffusion() const
  { return unitMass() * unitTime()/ unitVolume(); }

  /*! Helper function to generate the units of mutal diffusion
      outputted by the output plugins (see OPThermalConductivityE)*/
  inline double unitThermalCond() const
  { return unitk()/(unitLength() * unitTime()); }

  /*! Helper function to generate the units of ThermalDiffusion
      outputted by the output plugins (see OPThermalDiffusionE)*/
  inline double unitThermalDiffusion() const
  { return unitMass()/(unitLength() * unitTime()); }  

  /*! Helper function to generate the units of Viscosity
      outputted by the output plugins (see OPViscosityE)*/
  inline double unitViscosity() const
  { return 1.0 /(unitLength() * unitTime()); }
  
  /*! Helper function to generate the units of Viscosity
      outputted by the output plugins (see OPViscosityE)*/
  inline double unitPressure() const
  { return unitMass() / (unitLength()*unitTime()*unitTime()); }

  /*! \brief Used to rescale the system size after a system compression*/
  virtual void rescaleLength(double) = 0;

  /*! \brief Calculates the volume of the system*/
  double simVolume() const;

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Units&);
  
  virtual void operator<<(const XMLNode &) = 0;
  
  static Units* loadUnits(const XMLNode&, const DYNAMO::SimData*);
  
 protected:
  virtual void outputXML(xml::XmlStream &) const = 0;
};

#endif

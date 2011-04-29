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
#include "../../base/constants.hpp"
#include <magnet/units.hpp>
#include <cmath>

namespace xml { class XmlStream; }
namespace magnet { namespace xml { class Node; } }

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
class Units
{
 public:
  Units(double unitLength = 1, double unitTime = 1):
    _unitLength(unitLength),
    _unitTime(unitTime)
  {}

  /*! \brief Returns the simulation unit of time */
  double unitTime() const { return _unitTime; }

  /*! \brief Returns the simulation unit of length */
  double unitLength() const { return _unitLength; }

  /*! \brief Sets the unit length */
  void setUnitLength(double ul) { _unitLength = ul; }

  /*! \brief Sets the unit length */
  void setUnitTime(double ut) { _unitTime = ut; }

  /*! \brief Simulation unit of mass */
  double unitMass() const { return 1.0; }
  
  /*! \brief Boltzmanns constant */
  double unitk() const { return 1.0; }

  inline double getScaling(const magnet::units::Units& units) const
  {
    return std::pow(unitLength(), units.getUnitsPower(magnet::units::Units::L))
      * std::pow(unitTime(), units.getUnitsPower(magnet::units::Units::T))
      * std::pow(unitMass(), units.getUnitsPower(magnet::units::Units::M));
  }
  
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
  inline void rescaleLength(double r) { _unitLength *= r; }

  inline friend xml::XmlStream& operator<<(xml::XmlStream& XML, const Units& u)
  { u.outputXML(XML); return XML; }
  
  inline void operator<<(const magnet::xml::Node&) {}
  
 protected:
  inline void outputXML(xml::XmlStream &) const {}
  double _unitLength;
  double _unitTime;
};

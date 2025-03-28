/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
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
#include <cmath>
#include <magnet/units.hpp>

namespace xml {
class XmlStream;
}
namespace magnet {
namespace xml {
class Node;
}
} // namespace magnet

namespace dynamo {
/*! \brief The class used to convert in and out of simulation units.

  This class tracks the units of the simulation and provides helper
  functions to generate more complicated units from the elementary
  mass/length scales.

  The purpose of this class is to make it very easy to perform
  simulations in a "computationally easy" set of units, but have the
  simulator input and output in another set of units.

  Examples of this are:
  - If the system size can be rescaled to a 1x1x1 box, applying
  periodic boundary conditions becomes a simple rounding operation
  (this feature has been removed from DynamO, as it greatly
  complicates the code for very little speedup).
  - In replica-exchange simulations, the configurations may all be
  run at a reduced temperature so that the velocities do not have to
  be rescaled during exchanges (not implemented yet in DynamO).

  This class has a particular design specification in that it must
  be fully initialised on construction, so that other classes may
  directly start conversions on the loading of configurations.
 */
class Units {
public:
  Units(double unitLength = 1, double unitTime = 1)
      : _unitLength(unitLength), _unitTime(unitTime) {}

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

  inline double getScaling(const magnet::units::Units &units) const {
    return std::pow(unitLength(),
                    units.getUnitsPower(magnet::units::Units::L)) *
           std::pow(unitTime(), units.getUnitsPower(magnet::units::Units::T)) *
           std::pow(unitMass(), units.getUnitsPower(magnet::units::Units::M));
  }

  /*! Helper function to generate the unit of velocity*/
  inline double unitVelocity() const { return unitLength() / unitTime(); }
  inline double unitAcceleration() const {
    return unitLength() / (unitTime() * unitTime());
  }
  /*! Helper function to generate the unit of energy*/
  inline double unitEnergy() const {
    return unitMass() * unitVelocity() * unitVelocity();
  }
  /*! Helper function to generate the unit of area*/
  inline double unitArea() const { return unitLength() * unitLength(); }
  /*! Helper function to generate the unit of volume*/
  inline double unitVolume() const {
    return unitLength() * unitLength() * unitLength();
  }
  /*! Helper function to generate the unit of momentum*/
  inline double unitMomentum() const { return unitMass() * unitVelocity(); }
  inline double unitInertia() const { return unitArea() * unitMass(); }

  // Some dimensions of some properties
  /*! Helper function to generate the units of diffusion outputted by
    the output plugins (see OPMSD/OPMSDCorrelator)*/
  inline double unitDiffusion() const { return unitArea() / unitTime(); }

  /*! Helper function to generate the units of mutal diffusion
    outputted by the output plugins (see
    OPMutualDiffusionE/OPMutualDiffusionGK)*/
  inline double unitMutualDiffusion() const {
    return unitMass() * unitTime() / unitVolume();
  }

  /*! Helper function to generate the units of mutal diffusion
    outputted by the output plugins (see OPThermalConductivityE)*/
  inline double unitThermalCond() const {
    return unitk() / (unitLength() * unitTime());
  }

  /*! Helper function to generate the units of ThermalDiffusion
    outputted by the output plugins (see OPThermalDiffusionE)*/
  inline double unitThermalDiffusion() const {
    return unitMass() / (unitLength() * unitTime());
  }

  /*! Helper function to generate the units of Viscosity
    outputted by the output plugins (see OPViscosityE)*/
  inline double unitViscosity() const {
    return 1.0 / (unitLength() * unitTime());
  }

  /*! Helper function to generate the units of Viscosity
    outputted by the output plugins (see OPViscosityE)*/
  inline double unitPressure() const {
    return unitMass() / (unitLength() * unitTime() * unitTime());
  }

  /*! \brief Used to rescale the length scale after a system compression*/
  inline void rescaleLength(double r) { _unitLength *= r; }

  /*! \brief Used to rescale the time scale after a system compression.
   *
   * This rescaling is done in proportion to the length rescale, so
   * that the energy and velocity scales are unchanged.
   */
  inline void rescaleTime(double r) { _unitTime *= r; }

  inline friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                                   const Units &u) {
    u.outputXML(XML);
    return XML;
  }

  inline void operator<<(const magnet::xml::Node &) {}

protected:
  inline void outputXML(magnet::xml::XmlStream &) const {}
  double _unitLength;
  double _unitTime;
};
} // namespace dynamo

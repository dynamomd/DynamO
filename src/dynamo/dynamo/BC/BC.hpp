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
/*! \file BC.hpp
 *  Defines the BC class.
 */

#pragma once
#include <dynamo/base.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
namespace xml {
class Node;
class XmlStream;
} // namespace xml
} // namespace magnet

namespace dynamo {
class IntEvent;
class Simulation;
class Particle;

/*! \brief The base class for the Boundary Conditions of the simulation.

  This class has a couple of partial specialisations for CSqBC
  (square) and CRectBC (rectangular) periodic boundary conditions.
  These are utilised by the BCSquarePeriodic and BCRectangularPeriodic periodic
  boundary condition classes. There is the infinite system case BCNone.  More
  exotic conditions are the shearing BCSquareLeesEdwards and
  BCRectangularLeesEdwards (Lees-Edwards) boundary condition and one for
  studying confined systems in the x direction BCSquarePeriodicExceptX.
 */
class BoundaryCondition : public dynamo::SimBase_const {
public:
  /*! \brief Just the constructor

    \param SD The Simulation pointer.
    \param aName The name of the overall class.
   */
  BoundaryCondition(const dynamo::Simulation *const SD, const char *aName)
      : SimBase_const(SD, aName) {}

  virtual ~BoundaryCondition() {};

  /*! \brief This determines the minimum image length of a position vector

    This will turn the coordinates of a particle into the coordinates
    of the primary simulation image. For relative position vectors
    this will give the minimum image vector.

    \param pos The position vector to be altered.
   */
  virtual void applyBC(Vector &pos) const = 0;

  /*! \brief This determines the minimum image length of a position
    vector and the adjusted velocity vector.

    Exactly the same as the applyBC(Vector  &) function except if a
    velocity altering is required as part of the boundary condition
    then this is done too. This is used by boundary conditions such
    as BCSquareLeesEdwards.

    \param pos The position vector to affect.
    \param vel The corresponding velocity vector to affect.
   */
  virtual void applyBC(Vector &pos, Vector &vel) const = 0;

  /*! \brief A predictive boundary condition.

     This returns the rounding of the vector carried out as though it
     was performed dt in the future. Used in predicting cell
     transitions across the simulation boundaries.

    This is used by BC's like BCSquareLeesEdwards.

    \param pos The position vector to affect.
    \param dt The time difference to predict at.
   */
  virtual void applyBC(Vector &pos, const double &dt) const = 0;

  /*! \brief Stream the boundary conditions forward in time.*/
  virtual void update(const double &) {};

  /*! \brief Load the Boundary condition from an XML file. */
  virtual void operator<<(const magnet::xml::Node &) = 0;

  /*! \brief A helper for writing BoundaryCondition's to an XmlStream. */
  friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &,
                                            const BoundaryCondition &);

  /*! \brief The class loader for boundary conditions. */
  static shared_ptr<BoundaryCondition> getClass(const magnet::xml::Node &,
                                                dynamo::Simulation *);

  double getDistance(const Particle &p1, const Particle &p2);

protected:
  /*! \brief The XML output for a BoundaryCondition class*/
  virtual void outputXML(magnet::xml::XmlStream &XML) const = 0;
};
} // namespace dynamo

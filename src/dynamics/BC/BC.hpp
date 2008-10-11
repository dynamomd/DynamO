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
/*! \file BC.hpp
 *  Defines the BC class.
 */

#ifndef BC_H
#define BC_H

#include "../../datatypes/vector.hpp"

#include "../../base/is_base.hpp"

class XMLNode;
namespace xmlw
{
class XmlStream;
}

class CIntEvent;
namespace DYNAMO
{
  class SimData;
}

/*! \brief The base class for the Boundary Conditions of the simulation.
 *
 * This class has a couple of partial specialisations for CSqBC
 * (square) and CRectBC (rectangular) periodic boundary conditions.
 * These are utilised by the CSPBC and CRPBC periodic boundary
 * condition classes. There is the infinite system case CNullBC.  More
 * exotic conditions are the shearing CSLEBC and CRLEBC (Lees-Edwards)
 * boundary condition and one for studying confined systems in the x
 * direction CRNoXPBC.
 */
class CBC: public DYNAMO::SimBase_const
{
 public:
  /*! \brief Just the constructor
   *
   * \param SD The SimData pointer.
   * \param aName The name of the overall class.
   * \param aColor The colour of the class output.
   */
  CBC(const DYNAMO::SimData* const SD, const char *aName, 
      const char *aColor):
    SimBase_const(SD, aName, aColor)
  {}

  virtual ~CBC() {};

  /*! \brief This determines the minimum image length of a position vector
   *
   * This will turn the coordinates of a particle into the coordinates
   * of the primary simulation image. For relative position vectors
   * this will give the minimum image vector.
   *
   * \param pos The position vector to be altered.
   */
  virtual void setPBC(CVector<> & pos)const = 0;

  /*! \brief This determines the minimum image length of a position
   * vector and the adjusted velocity vector.
   *
   * Exactly the same as the setPBC(CVector<> &) function except if a
   * velocity altering is required as part of the boundary condition
   * then this is done too. This is used by boundary conditions such
   * as CSLEBC.
   *
   * \param pos The position vector to affect.
   * \param vel The corresponding velocity vector to affect.
   */
  virtual void setPBC(CVector<> & pos, CVector<> & vel) const = 0;

  /*! \brief A predictive boundary condition.
   *
   *  This returns the rounding of the vector carried out as though it
   *  was performed dt in the future. Used in predicting cell
   *  transitions across the simulation boundaries.
   *
   * This is used by BC's like CSLEBC.
   *
   * \param pos The position vector to affect.
   * \param dt The time difference to predict at.
   */
  virtual void setPBC(CVector<> &pos, Iflt dt) const 
    { setPBC(pos); }

  /*! \brief Stream the boundary conditions forward in time.*/
  virtual void update(const Iflt&) {};

  /*! \brief Load the Boundary condition from an XML file. */
  virtual void operator<<(const XMLNode &) = 0;

  /*! \brief A polymorphic class copy helper. */
  virtual CBC* Clone () const = 0;

  /*! \brief A helper for writing CBC's to an XmlStream. */
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CBC&);

  /*! \brief The class loader for boundary conditions. */
  static CBC* loadClass(const XMLNode& ,DYNAMO::SimData*);

 protected:
  /*! \brief The XML output for a CBC class*/
  virtual void outputXML(xmlw::XmlStream &XML) const = 0;

  /*! \brief The actual positional rounding to carry out on a position
   * vector.
   */
  virtual void rounding(CVector<> &) const = 0;
};

#endif

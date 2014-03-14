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
#include <dynamo/base.hpp>
#include <dynamo/ranges/IDRange.hpp>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

namespace dynamo {
  class IntEvent;
  class NEventData;
  class GlobalEvent;

  /*! \brief Base class for Non-\ref Local single-particle events.
   *
   * A global event is a single particle event which cannot be optimized
   * by using a neighbour list. In fact, neighbour lists are Global
   * event types and have a specialization of the Global interface (\ref
   * GNeighbourList).
   */
  class Global: public dynamo::SimBase
  {
  public:
    /*! \brief Constructor.
     *
     * \param sim A pointer to the root of the simulation data.
     * \param name The name of the class (for formatted output).
     * \param range The range of particles for which this interaction is
     * valid (the default value of NULL indicates all particles are valid).
     */
    Global(dynamo::Simulation* sim, std::string name, IDRange* range = NULL);
  
    /*! \brief Returns true if the Global applies to the passed
     * particle.
     */
    bool isInteraction(const Particle&) const;

    /*! \brief Returns the next calculated event for the passed
     * particle.
     */
    virtual GlobalEvent getEvent(const Particle &) const = 0;

    /*! \brief Executes the event for a particle.
     * 
     * \param p The particle which is about to undergo an interaction.
     * \param dt The time the scheduler thinks this particles Global
     * event will occur in.
     */
    virtual void runEvent(Particle& p, const double dt) = 0;

    /*! \brief Initializes the Global event.
     */
    virtual void initialise(size_t nID)  { ID=nID; }

    /*! \brief Helper function for saving an XML representation of this
     * class.
     */
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Global&);

    /*! \brief Constructs a derived Global class according to the passed
     * XML Node.
     */
    static shared_ptr<Global> getClass(const magnet::xml::Node&, dynamo::Simulation*);

    /*! \brief Loads a derived class from a passed XML Node.
     */
    virtual void operator<<(const magnet::xml::Node&) = 0;

    /*! \brief Sets the name by which this Global is referred to.
     */
    void setName(const std::string& tmp) { globName = tmp; }

    /*! \brief Returns the name by which this Global is referred to.
     */
    const std::string& getName() const { return globName; }

    /*! \brief Returns the unique ID number of this Global.
     */
    inline const size_t& getID() const { return ID; }
  
  protected:
    /*! \brief Writes out an XML representation of the Global
     */
    virtual void outputXML(magnet::xml::XmlStream&) const = 0;

    shared_ptr<IDRange> range;  
    std::string globName;
    size_t ID;
  };
}


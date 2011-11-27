/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/base.hpp>
#include <dynamo/dynamics/ranges/2range.hpp>
#include <string>

namespace magnet { namespace xml { class Node; class XmlStream; } }

namespace dynamo {
  class Range;
  class PairEventData;
  class IntEvent;
  class Species;

  /*! \brief This class is the base interface for Interation classes.
  
   Interaction's are events that describe the Interaction between two
   particles. These classes are responsible for: 

   - Storing the values used in calculating the interactions (e.g.,
     the interaction diameter).

   - Storing the "state" of the interaction, to ensure only correct
     dynamics occur (e.g., a square well particle must capture a
     particle before it can be released or it can hit the inner
     core). This state storing often uses one of the ICapture classes
     (ISingleCapture or IMultiCapture).

   - Performing high level calculations or optimizations for the
     interactions (e.g., for hard lines (ILines), we use a bounding
     sphere before testing for the expensive line-line collision, the
     bounding sphere test is organised in here).

   \warning You must only perform high level calculations here. All
   actual collision testing should use the "primative" functions
   defined in the Liouvillean class. This allows an interaction to be
   easily ported to alternative dynamics (like compression or
   gravity).
  */
  class Interaction: public dynamo::SimBase
  {
  public:
    Interaction(dynamo::SimData*, C2Range*);
  
    virtual ~Interaction() {}

    virtual void initialise(size_t) = 0;

    //! Calculate if an event is to occur between two particles.
    virtual IntEvent getEvent(const Particle &, 
			      const Particle &) const = 0;

    //! Run the dynamics of an event that is occuring now.
    virtual void runEvent(const Particle&, const Particle&, const IntEvent&) const = 0;

    //! Return the maximum distance at which two particles may interact using this Interaction.
    //!
    //! This value is used in GNeighbourList's to make sure a certain GNeighbourList is suitable for detecting possible Interaction partner particles.
    virtual double maxIntDist() const = 0;  

    //! Returns the internal energy "stored" in this interaction.
    virtual double getInternalEnergy() const = 0; 

    //! Returns the internal energy "stored" in this interaction by the
    //! two passed particles.
    virtual double getInternalEnergy(const Particle&, const Particle&) const { return 0; } 

    //! Returns the excluded volume of a certain particle.
    virtual double getExcludedVolume(size_t) const = 0;

    //! Loads the parameters of the Interaction from an XML node in the configuration file.
    virtual void operator<<(const magnet::xml::Node&);
  
    //! A helper function that calls Interaction::outputXML to write out the parameters of this interaction to a config file.
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Interaction&);
 
    //! This static function will instantiate a new interaction of the correct type specified by the xml node passed. 
    //!
    //! This is the birth point for all Interactions loaded from a configuration file.
    static shared_ptr<Interaction> getClass(const magnet::xml::Node&, dynamo::SimData*);

    //! Tests if this interaction is meant to be used between the two passed Particle -s.
    bool isInteraction(const Particle &p1, const Particle &p2) const
    { return range->isInRange(p1,p2); }
  
    //! Tests if this interaction may have been used for the passed interaction event (IntEvent).
    bool isInteraction(const IntEvent &) const;

    //! Tests if this interaction is suitable to describe the basic properties of an entire species.
    bool isInteraction(const Species &) const;

    //! Returns the "name" of the interaction used in name-based look-ups.
    inline const std::string& getName() const { return intName; }

    //! Returns the C2Range describing the pairs of particles this
    //! Interaction can generate events for.
    shared_ptr<C2Range>& getRange();

    //! Returns the C2Range describing the pairs of particles this
    //! Interaction can generate events for.
    const shared_ptr<C2Range>& getRange() const;

    //! Test if an invalid state has occurred between the two passed particles
    virtual void checkOverlaps(const Particle&, const Particle&) const = 0;

    //! Return the ID number of the Interaction. Used for fast look-ups,
    //! once a name-based look up has been completed.
    inline const size_t& getID() const { return ID; }

  protected:
    //! This constructor is only to be used when using virtual
    //! inheritance, the bottom derived class must explicitly call the
    //! other Interaction ctor.
    Interaction() { M_throw() << "Default constructor called!"; }

    //! Write out an XML tag that describes this Interaction and stores its Property -s.
    virtual void outputXML(magnet::xml::XmlStream&) const = 0;

    shared_ptr<C2Range> range;

    std::string intName;
    size_t ID;
  };
}


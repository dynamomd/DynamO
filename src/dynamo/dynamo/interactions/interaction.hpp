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

#include <array>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/base.hpp>
#include <dynamo/eventtypes.hpp>
#include <dynamo/ranges/IDPairRange.hpp>
#include <limits>
#include <string>

namespace magnet {
namespace xml {
class Node;
class XmlStream;
} // namespace xml
} // namespace magnet

namespace dynamo {
class IDRange;
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
 defined in the Dynamics class. This allows an interaction to be
 easily ported to alternative dynamics (like compression or
 gravity).
*/
class Interaction : public dynamo::SimBase {
public:
  Interaction(dynamo::Simulation *, IDPairRange *);

  virtual ~Interaction() {}

  virtual void initialise(size_t id) {
    ID = id;
    setOutputPrefix(std::string("Interaction[") + intName + std::string("]"));
  }

  /*! \brief Calculate if and when an event is to occur between two
      particles.
   */
  virtual Event getEvent(const Particle &, const Particle &) const = 0;

  /*! \brief Run the dynamics of an event which is occuring now.
   */
  virtual PairEventData runEvent(Particle &, Particle &, Event) = 0;

  /*! \brief Return the maximum distance at which two particles may interact
    using this Interaction.

    This value is used in GNeighbourList's to make sure a certain
    GNeighbourList is suitable for detecting possible Interaction
    partner particles.
  */
  virtual double maxIntDist() const = 0;

  /*! \brief Returns the internal energy "stored" in this interaction.
   */
  virtual double getInternalEnergy() const { return 0; }

  /*! \brief Returns the internal energy "stored" in this
   interaction by the two passed particles.
  */
  virtual double getInternalEnergy(const Particle &, const Particle &) const {
    return 0;
  }

  /*! \brief Returns the excluded volume of a certain particle.
   */
  virtual double getExcludedVolume(size_t) const = 0;

  /*! \brief Loads the parameters of the Interaction from an XML
      node in the configuration file.
   */
  virtual void operator<<(const magnet::xml::Node &);

  /*! \brief A helper function that calls Interaction::outputXML to
      write out the parameters of this interaction to a config file.
   */
  friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &,
                                            const Interaction &);

  /*! \brief This static function will instantiate a new interaction
      of the correct type specified by the xml node passed.

    This is the birth point for all Interactions loaded from a
    configuration file.
  */
  static shared_ptr<Interaction> getClass(const magnet::xml::Node &,
                                          dynamo::Simulation *);

  /*! \brief Tests if this interaction is meant to be used between
      the two passed Particle -s.
   */
  bool isInteraction(const Particle &p1, const Particle &p2) const {
    return range->isInRange(p1, p2);
  }

  /*! \brief Tests if this interaction may have been used for the
      passed interaction event (IntEvent).
   */
  bool isInteraction(const Event &) const;

  /*! \brief Returns the "name" of the interaction used in
      name-based look-ups.
   */
  inline const std::string &getName() const { return intName; }

  /*! \brief Returns the IDPairRange describing the pairs of particles
   this Interaction can generate events for.
  */
  shared_ptr<IDPairRange> &getRange();

  /*! \brief Returns the IDPairRange describing the pairs of particles
      this Interaction can generate events for.
   */
  const shared_ptr<IDPairRange> &getRange() const;

  /*! \brief Test if an invalid state has occurred between the two
      passed particles.

      \param p1 The first particle of the pair to test.
      \param p2 The second particle of the pair to test.
      \param textoutput If true, there will be a text description of
      any detected invalid states printed to cerr.
      \return true if the pair is in an invalid state.
  */
  virtual bool validateState(const Particle &p1, const Particle &p2,
                             bool textoutput = true) const = 0;

  /*! \brief Test if the internal state of the Interaction is valid.

      \param textoutput If true, there will be a text description of
      any detected invalid states printed to cerr.
      \param max_reports The maximum number of invalid state reports to write to
     cerr.
      \return The total count of the invalid internal states.
  */
  virtual size_t
  validateState(bool textoutput = true,
                size_t max_reports = std::numeric_limits<size_t>::max()) const {
    return 0;
  }

  /*! \brief Return the ID number of the Interaction. Used for fast
   look-ups, once a name-based look up has been completed.
  */
  inline const size_t &getID() const { return ID; }

  virtual void outputData(magnet::xml::XmlStream &) const {}

  enum GLYPH_TYPE {
    SPHERE_GLYPH = 0,
    ARROW_GLYPH = 1,
    CYLINDER_GLYPH = 2,
    ROD_GLYPH = 3,
    LINE_GLYPH = 4,
    CUBE_GLYPH = 5,
    DUMBBELL_GLYPH = 6
  };

  /*! Returns the diameter of the one of the spheres used to represent
    a particle.

    \param ID The id of the particle to fetch the diameter for.

    \param subID The index of one of the spheres used to render the
    particle. This must be in the range [0, spheresPerParticle()-1].
  */
  virtual std::array<double, 4> getGlyphSize(size_t ID) const {
    M_throw() << "Cannot determine the glyph size of particle " << ID
              << " using the Interaction \"" << getName() << "\"."
              << "To visualise this system, please make sure the "
                 "self-Interactions of each particle are drawable.";
  }

  /*! Returns the location of the one of the spheres used to
    represent a particle.

    \param ID The id of the particle to fetch the location for.

    \param subID The index of one of the spheres used to render the
    particle. This must be in the range [0, spheresPerParticle()-1].
  */
  virtual GLYPH_TYPE getDefaultGlyphType() const { return SPHERE_GLYPH; }

protected:
  /*! \brief This constructor is only to be used when using virtual
   inheritance, the bottom derived class must explicitly call the
   other Interaction constructor.
  */
  Interaction() { M_throw() << "Default constructor called!"; }

  /*! \brief Write out an XML tag that describes this Interaction
      and stores its Property -s.
   */
  virtual void outputXML(magnet::xml::XmlStream &) const = 0;

  shared_ptr<IDPairRange> range;

  std::string intName;
  size_t ID;
};
} // namespace dynamo

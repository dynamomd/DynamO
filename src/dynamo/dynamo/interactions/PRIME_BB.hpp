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

#include <dynamo/interactions/interaction.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/PRIME.hpp>
#include <dynamo/interactions/captures.hpp>

namespace dynamo {
  class IPRIME_BB: public ICapture
  {
  public:
    IPRIME_BB(const magnet::xml::Node&, dynamo::Simulation*);

    virtual std::array<double, 4> getGlyphSize(size_t ID) const;

    virtual double getExcludedVolume(size_t) const;

    virtual double maxIntDist() const;

    virtual void initialise(size_t);

    virtual size_t captureTest(const Particle&, const Particle&) const;

    virtual IntEvent getEvent(const Particle&, const Particle&) const;

    virtual PairEventData runEvent(Particle&, Particle&, const IntEvent&);

    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual double getInternalEnergy(const Particle&, const Particle&) const;

    virtual double getInternalEnergy() const;

    virtual void operator<<(const magnet::xml::Node&);

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

  protected:
    /*! \brief Returns the type of the bead on the backbone.
     */
    TPRIME::BeadData getBeadData(const size_t particleID) const {
      return _topology->getBeadInfo(particleID);
    }

    /*! \brief Calculates the interaction parameters for the passed pair.

      \return This pair has the interaction diameter as the first
      value and whether this diameter is a bond as the second value.
     */
    std::tuple<double, double, double> getInteractionParameters(const size_t pID1, const size_t pID2) const;

    std::shared_ptr<TPRIME> _topology;

    double _PRIME_HB_strength;
  };
}

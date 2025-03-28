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
#include <magnet/math/matrix.hpp>

namespace dynamo {
class IParallelCubes : public Interaction {
public:
  template <class T1, class T2>
  IParallelCubes(dynamo::Simulation *tmp, T1 d, T2 e, IDPairRange *nR,
                 std::string name)
      : Interaction(tmp, nR),
        _diameter(Sim->_properties.getProperty(d, Property::Units::Length())),
        _e(Sim->_properties.getProperty(e, Property::Units::Dimensionless())) {
    intName = name;
  }

  IParallelCubes(const magnet::xml::Node &, dynamo::Simulation *);

  virtual std::array<double, 4> getGlyphSize(size_t ID) const;
  virtual GLYPH_TYPE getDefaultGlyphType() const { return CUBE_GLYPH; }

  void operator<<(const magnet::xml::Node &);

  virtual double maxIntDist() const;

  virtual double getExcludedVolume(size_t) const;

  virtual Event getEvent(const Particle &, const Particle &) const;

  virtual PairEventData runEvent(Particle &, Particle &, Event);

  virtual void outputXML(magnet::xml::XmlStream &) const;

  virtual bool validateState(const Particle &p1, const Particle &p2,
                             bool textoutput = true) const;

protected:
  shared_ptr<Property> _diameter;
  shared_ptr<Property> _e;
};
} // namespace dynamo

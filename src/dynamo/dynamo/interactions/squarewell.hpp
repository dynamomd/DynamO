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

#include <dynamo/interactions/captures.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
class ISquareWell : public ICapture {
public:
  template <class T1, class T2, class T3, class T4>
  ISquareWell(dynamo::Simulation *tmp, T1 d, T2 l, T3 wd, T4 e, IDPairRange *nR,
              std::string name)
      : ICapture(tmp, nR),
        _diameter(Sim->_properties.getProperty(d, Property::Units::Length())),
        _lambda(
            Sim->_properties.getProperty(l, Property::Units::Dimensionless())),
        _wellDepth(Sim->_properties.getProperty(wd, Property::Units::Energy())),
        _e(Sim->_properties.getProperty(e, Property::Units::Dimensionless())) {
    intName = name;
  }

  ISquareWell(const magnet::xml::Node &, dynamo::Simulation *);

  void operator<<(const magnet::xml::Node &);

  virtual std::array<double, 4> getGlyphSize(size_t ID) const;

  virtual double getExcludedVolume(size_t) const;

  virtual double maxIntDist() const;

  virtual size_t captureTest(const Particle &, const Particle &) const;

  virtual void initialise(size_t);

  virtual Event getEvent(const Particle &, const Particle &) const;

  virtual PairEventData runEvent(Particle &, Particle &, Event);

  virtual void outputXML(magnet::xml::XmlStream &) const;

  virtual double getInternalEnergy(const Particle &, const Particle &) const;
  virtual double getInternalEnergy() const;

  using ICapture::validateState;
  virtual bool validateState(const Particle &p1, const Particle &p2,
                             bool textoutput = true) const;

protected:
  ISquareWell(dynamo::Simulation *tmp, IDPairRange *nR) : ICapture(tmp, nR) {}

  shared_ptr<Property> _diameter;
  shared_ptr<Property> _lambda;
  shared_ptr<Property> _wellDepth;
  shared_ptr<Property> _e;
};
} // namespace dynamo

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
class IDumbbells : public ICapture {
public:
  template <class T1>
  IDumbbells(dynamo::Simulation *tmp, T1 e, IDPairRange *nR, std::string name)
      : ICapture(tmp, nR),
        _e(Sim->_properties.getProperty(e, Property::Units::Dimensionless())),
        _unusedDimension(std::numeric_limits<size_t>::max()) {
    intName = name;
  }

  virtual std::array<double, 4> getGlyphSize(size_t ID) const;

  virtual GLYPH_TYPE getDefaultGlyphType() const { return DUMBBELL_GLYPH; }

  virtual double getExcludedVolume(size_t ID) const;

  IDumbbells(const magnet::xml::Node &, dynamo::Simulation *);

  void operator<<(const magnet::xml::Node &);

  virtual void initialise(size_t);

  virtual double maxIntDist() const;

  double maxIntDist(size_t p1, size_t p2) const;

  virtual Event getEvent(const Particle &, const Particle &) const;

  virtual PairEventData runEvent(Particle &, Particle &, Event);

  virtual void outputXML(magnet::xml::XmlStream &) const;

  virtual size_t captureTest(const Particle &, const Particle &) const;

  using ICapture::validateState;
  virtual bool validateState(const Particle &p1, const Particle &p2,
                             bool textoutput = true) const;

  void setUnusedDimension(size_t v) { _unusedDimension = v; }

  template <class T1> void addSphere(const Vector &offset, const T1 &diam) {
    _compositeData.push_back(Sphere(
        offset, Sim->_properties.getProperty(diam, Property::Units::Length()),
        Sim->_properties.getProperty(0, Property::Units::Energy()),
        Sim->_properties.getProperty(0, Property::Units::Dimensionless())));
  }

protected:
  struct Sphere {
    Sphere(Vector offset, shared_ptr<Property> diam,
           shared_ptr<Property> welldepth, shared_ptr<Property> lambda)
        : _offset(offset), _diam(diam), _welldepth(welldepth), _lambda(lambda) {
    }
    Vector _offset;
    shared_ptr<Property> _diam;
    shared_ptr<Property> _welldepth;
    shared_ptr<Property> _lambda;
  };

  std::vector<Sphere> _compositeData;
  shared_ptr<Property> _e;
  size_t _unusedDimension;
};
} // namespace dynamo

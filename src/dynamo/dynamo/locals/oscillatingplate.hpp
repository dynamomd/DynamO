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

#include <dynamo/coilRenderObj.hpp>
#include <dynamo/locals/local.hpp>

namespace coil {
class RSurface;
}
namespace dynamo {
class LOscillatingPlate : public Local, public CoilRenderObj {
public:
  LOscillatingPlate(const magnet::xml::Node &, dynamo::Simulation *);
  LOscillatingPlate(dynamo::Simulation *, Vector, Vector, double, double,
                    double, double, double, std::string, IDRange *,
                    double timeshift = 0, bool nstrongPlate = false);

  virtual ~LOscillatingPlate() {}

  virtual Event getEvent(const Particle &) const;

  virtual ParticleEventData runEvent(Particle &, const Event &) const;

  virtual void operator<<(const magnet::xml::Node &);

  Vector getPosition() const;

  Vector getVelocity() const;

  double getPlateEnergy() const;

  const Vector &getCentre() const { return rw0; }

  virtual bool validateState(const Particle &part,
                             bool textoutput = true) const {
    return false;
  }

#ifdef DYNAMO_visualizer
  virtual shared_ptr<coil::RenderObj> getCoilRenderObj() const;
  virtual void updateRenderData() const;
#endif

protected:
#ifdef DYNAMO_visualizer
  mutable shared_ptr<coil::RSurface> _renderObj;
#endif

  virtual void outputXML(magnet::xml::XmlStream &) const;

  bool strongPlate;
  Vector rw0;
  Vector nhat;
  double omega0;
  double sigma;
  double e;
  mutable double delta;
  double mass;
  mutable double timeshift;
  mutable size_t lastID;
  mutable long double lastsystemTime;
};
} // namespace dynamo

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

#include <dynamo/species/inertia.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
class SpFixedCollider : public SpInertia {
public:
  SpFixedCollider(dynamo::Simulation *sim, IDRange *r, std::string nName,
                  unsigned int ID)
      : SpInertia(sim, r, std::numeric_limits<float>::infinity(), nName, ID) {
    setOutputPrefix("SpFixedCollider");
  }

  SpFixedCollider(const magnet::xml::Node &XML, dynamo::Simulation *nSim,
                  unsigned int nID)
      : SpInertia(nSim, NULL, std::numeric_limits<float>::infinity(), "", nID) {
    setOutputPrefix("SpFixedCollider");
    operator<<(XML);
  }

  virtual void initialise();

  virtual void operator<<(const magnet::xml::Node &XML);

  virtual double getScalarMomentOfInertia(size_t ID) const {
    return std::numeric_limits<float>::infinity();
  }

  virtual double getParticleKineticEnergy(size_t ID) const { return 0; }

  virtual double getDOF() const { return 0; }

protected:
  virtual void outputXML(magnet::xml::XmlStream &XML) const;
};
} // namespace dynamo

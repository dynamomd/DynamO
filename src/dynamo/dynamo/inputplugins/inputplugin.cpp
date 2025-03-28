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

#include <dynamo/BC/include.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/inputplugins/inputplugin.hpp>
#include <dynamo/interactions/include.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/topology/include.hpp>

namespace dynamo {
InputPlugin::InputPlugin(dynamo::Simulation *tmp, const char *aName)
    : SimBase(tmp, aName) {}

void InputPlugin::rescaleVels(double val) {
  dout << "WARNING Rescaling kT to " << val << std::endl;

  double currentkT(Sim->dynamics->getkT() / Sim->units.unitEnergy());

  dout << "Current kT " << currentkT << std::endl;

  Sim->dynamics->rescaleSystemKineticEnergy(val / currentkT);
}

void InputPlugin::setCOMVelocity(const Vector vel) {
  dout << "Setting COM Velocity" << std::endl;

  if (Sim->N() <= 1)
    derr << "Refusing to set momentum for a " << Sim->N() << " particle system"
         << std::endl;
  else
    Sim->setCOMVelocity(vel);
}

void InputPlugin::zeroMomentum() {
  dout << "Zeroing Momentum" << std::endl;

  if (Sim->N() <= 1)
    derr << "Refusing to zero momentum for a " << Sim->N() << " particle system"
         << std::endl;
  else
    Sim->setCOMVelocity();
}

void InputPlugin::zeroCentreOfMass() {
  dout << "Zeroing Centre of Mass" << std::endl;

  Vector com{0, 0, 0};
  double totmass = 0.0;
  for (const Particle &part : Sim->particles) {
    totmass += Sim->species[part]->getMass(part);
    com += part.getPosition() * Sim->species[part]->getMass(part);
  }
  com /= totmass;

  for (Particle &part : Sim->particles)
    part.getPosition() -= com;
}

void InputPlugin::mirrorDirection(unsigned int iDim) {
  for (Particle &part : Sim->particles) {
    part.getVelocity()[iDim] *= -1.0;
    part.getPosition()[iDim] *= -1.0;
  }
}

void InputPlugin::zeroVelComp(size_t iDim) {
  dout << "Zeroing the " << iDim << " dimension velocities" << std::endl;
  for (Particle &part : Sim->particles)
    part.getVelocity()[iDim] = 0.0;
}
} // namespace dynamo

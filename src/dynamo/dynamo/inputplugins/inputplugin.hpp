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
#include <magnet/math/vector.hpp>

namespace dynamo {
class InputPlugin : public dynamo::SimBase {
public:
  InputPlugin(dynamo::Simulation *, const char *aName);

  virtual ~InputPlugin() {};

  virtual void initialise() {};

  // Rescaling system
  void rescaleVels(double val = 1.0);

  void zeroMomentum();

  void setCOMVelocity(const Vector);

  void zeroCentreOfMass();

  void setPackFrac(double);

  void mirrorDirection(unsigned int);

  void zeroVelComp(size_t);

protected:
};
} // namespace dynamo

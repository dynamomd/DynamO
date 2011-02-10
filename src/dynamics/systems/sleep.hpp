/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include "system.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include <magnet/cloneptr.hpp>
#include <map>

class SSleep: public System
{
public:
  SSleep(const XMLNode& XML, DYNAMO::SimData*);

  SSleep(DYNAMO::SimData*, std::string, CRange*, double);
  
  virtual System* Clone() const { return new SSleep(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xml::XmlStream&) const;

  void particlesUpdated(const NEventData&);

  void recalculateTime();

  bool sleepCondition(const Particle& part, const Vector& g, const Vector& vel = Vector(0,0,0));

  magnet::ClonePtr<CRange> _range;
  double _sleepDistance;
  double _sleepTime;
  double _sleepVelocity;

  mutable std::map<size_t, Vector> stateChange;

  std::vector<std::pair<Vector, long double> > _lastData;
};

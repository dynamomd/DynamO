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

#include <dynamo/2particleEventData.hpp>
#include <list>

namespace dynamo {
  class NEventData
  {
  public:
    NEventData() {};
    NEventData(const ParticleEventData& a) { L1partChanges.push_back(a); }
    NEventData(const PairEventData& a) { L2partChanges.push_back(a); }

    NEventData&  operator+=(const ParticleEventData& p) { L1partChanges.push_back(p); return *this; }
    NEventData&  operator+=(const PairEventData& p) { L2partChanges.push_back(p); return *this; }

    std::vector<ParticleEventData> L1partChanges;
    std::vector<PairEventData> L2partChanges;
  };
}

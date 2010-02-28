/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CNParticleData_HPP
#define CNParticleData_HPP

#include "2particleEventData.hpp"
#include <list>

class CNParticleData
{
public:
  CNParticleData() {};
  CNParticleData(const C1ParticleData& a) { L1partChanges.push_back(a); }
  CNParticleData(const C2ParticleData& a) { L2partChanges.push_back(a); }

  std::list<C1ParticleData> L1partChanges;
  std::list<C2ParticleData> L2partChanges;
};

#endif

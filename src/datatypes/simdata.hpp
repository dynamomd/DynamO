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

#ifndef SSimBase_H
#define SSimBase_H

#include <vector>

struct SSimBase
{
  SSimBase(  const std::vector<CParticle> &particleList2,
	     const Dynamics * const dynamics2,
	     const long &nColl2,
	     const long &maxColl2,
	     const long &nPrint2,
	     const Iflt &sysTime2,
	     const std::ostringstream &history2
	     ):
    particleList(particleList2), dynamics(dynamics2),
       nColl(nColl2), maxColl(maxColl2), nPrint(nPrint2),
       sysTime(sysTime2), history(history2) {};
    
  const std::vector<CParticle> &particleList;  
  const Dynamics * const dynamics;
  const long &nColl;
  const long &maxColl;
  const long &nPrint;
  const Iflt &sysTime;
  const std::ostringstream &history;

};

#endif

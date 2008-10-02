/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef SimImage_H
#define SimImage_H

#include <list>
#include <vector>
#include "../datatypes/pluginpointer.hpp"
#include "../dynamics/dynamics.hpp"
class CDynamics;
class COutputPlugin;
class CParticle;

class CSimImage
{
 public:
  CSimImage(Iflt, long, std::vector<CParticle>,
	    std::vector<smrtPlugPtr<COutputPlugin> >, const CDynamics&);

  CSimImage(const CSimImage&);
  
  Iflt sysTime;
  unsigned long nColl;
  std::vector<CParticle> particleList;
  std::vector<smrtPlugPtr<COutputPlugin> > COPlugins; 
  CDynamics dynamics;
};

#endif

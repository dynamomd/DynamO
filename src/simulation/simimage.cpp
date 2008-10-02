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

#include "simimage.hpp"
#include "../dynamics/dynamics.hpp"
#include "particle.hpp"
#include "../outputplugins/outputplugin.hpp"

CSimImage::CSimImage(Iflt t, long c, std::vector<CParticle> pList, 
		     std::vector<smrtPlugPtr<COutputPlugin> > COPs,const  CDynamics &dyn) 
  :sysTime(t), nColl(c), particleList(pList), COPlugins(COPs),
   dynamics(dyn)
{}

CSimImage::CSimImage(const CSimImage &SI)
  :sysTime(SI.sysTime), nColl(SI.nColl), 
   particleList(SI.particleList), COPlugins(SI.COPlugins),
   dynamics(SI.dynamics)
{}

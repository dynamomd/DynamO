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

#include "velprof.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"

COPVelProfile::COPVelProfile(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"VelProfile"),
  samplesTaken(0),
  binWidth(0.5)
{
  if (NDIM != 3)
    D_throw() << "Terrible plugin for 3 dims only";
}

void 
COPVelProfile::initialise()
{    
  binWidth *= Sim->Dynamics.units().unitLength();
  
  vx.resize(static_cast<size_t>(Sim->aspectRatio[1]/binWidth)+1,
	    std::vector<std::pair<size_t, Iflt> >
	    (static_cast<size_t>(Sim->aspectRatio[2]/binWidth)+1,
	     std::pair<size_t, Iflt>(0, 0)));

}

void 
COPVelProfile::ticker()
{
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      Vector  pos(part.getPosition());
      Vector  vel(part.getVelocity());

      Sim->Dynamics.BCs().applyBC(pos, vel);

      pos += Sim->aspectRatio / 2.0;

      std::pair<size_t, Iflt>& locpair(vx[static_cast<size_t>(pos[1]/binWidth)]
				       [static_cast<size_t>(pos[2]/binWidth)]);
      ++locpair.first;
      locpair.second += vel[0];
      
    }

}

void 
COPVelProfile::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("VelProfile")
      << xmlw::chardata();
  
  size_t nybins = static_cast<size_t>(Sim->aspectRatio[1]/binWidth)+1;
  size_t nzbins = static_cast<size_t>(Sim->aspectRatio[2]/binWidth)+1;

  for (size_t y = 0; y < nybins; ++y)
    {
      for (size_t z = 0; z < nzbins; ++z)
	XML << y * binWidth / Sim->Dynamics.units().unitLength()
	    << " " << z * binWidth / Sim->Dynamics.units().unitLength()
	    << " "
	    << ((vx[y][z].first) ? (vx[y][z].second / vx[y][z].first) 
		/ Sim->Dynamics.units().unitVelocity() : 0)
	    << "\n";

      XML << "\n";
    }
      
  XML << xmlw::endtag("VelProfile");
}

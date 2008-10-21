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

#include "vel_dist.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

COPVelDist::COPVelDist(const DYNAMO::SimData* tmp, 
		       const XMLNode& XML):
  COPTicker(tmp,"VelDist"),
  binWidth(0.01)
{ operator<<(XML); }

void 
COPVelDist::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("binWidth"))
      binWidth = boost::lexical_cast<Iflt>(XML.getAttribute("binWidth"));
      }
  catch (std::exception& excep)
    {
      D_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void 
COPVelDist::initialise()
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    data[iDim].resize(Sim->Dynamics.getSpecies().size(), 
		      C1DHistogram(Sim->Dynamics.units().unitVelocity() 
				   * binWidth));
}

void 
COPVelDist::ticker()
{
  Sim->Dynamics.Liouvillean().updateAllParticles();
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *sp.getRange())
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if (iDim == 0)
	data[iDim][sp.getID()]
	  .addVal(Sim->vParticleList[ID].getVelocity()[0] 
		  - Sim->vParticleList[ID].getPosition()[1]);
      else
	data[iDim][sp.getID()]
	  .addVal(Sim->vParticleList[ID].getVelocity()[iDim]);
}

void
COPVelDist::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("VelDist");
  
  for (size_t id = 0; id < Sim->Dynamics.getSpecies().size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->Dynamics.getSpecies()[id].getName();
     
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	{
	  XML << xmlw::tag("Dimension")
	      << xmlw::attr("val")
	      << iDim;
	  
	  data[iDim][id].outputHistogram
	    (XML, 1.0 / Sim->Dynamics.units().unitVelocity());
	  
	  XML << xmlw::endtag("Dimension");
	}
      
      XML << xmlw::endtag("Species");
    }
  
  XML << xmlw::endtag("VelDist");
}

/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

OPVelDist::OPVelDist(const DYNAMO::SimData* tmp, 
		       const XMLNode& XML):
  OPTicker(tmp,"VelDist"),
  binWidth(0.01)
{ operator<<(XML); }

void 
OPVelDist::operator<<(const XMLNode& XML)
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
OPVelDist::initialise()
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    data[iDim].resize(Sim->dynamics.getSpecies().size(), 
		      C1DHistogram(Sim->dynamics.units().unitVelocity() 
				   * binWidth));
}

void 
OPVelDist::ticker()
{
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *sp->getRange())
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      data[iDim][sp->getID()]
	.addVal(Sim->vParticleList[ID].getVelocity()[iDim]);
}

void
OPVelDist::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("VelDist");
  
  for (size_t id = 0; id < Sim->dynamics.getSpecies().size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->dynamics.getSpecies()[id]->getName();
     
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	{
	  XML << xmlw::tag("Dimension")
	      << xmlw::attr("val")
	      << iDim;
	  
	  data[iDim][id].outputHistogram
	    (XML, 1.0 / Sim->dynamics.units().unitVelocity());
	  
	  XML << xmlw::endtag("Dimension");
	}
      
      XML << xmlw::endtag("Species");
    }
  
  XML << xmlw::endtag("VelDist");
}

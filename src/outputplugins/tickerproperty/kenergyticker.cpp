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

#include "kenergyticker.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

COPKEnergyTicker::COPKEnergyTicker(const DYNAMO::SimData* tmp, 
		       const XMLNode& XML):
  COPTicker(tmp,"KEnergyTicker"),
  count(0)
{ operator<<(XML); }

void 
COPKEnergyTicker::operator<<(const XMLNode& XML)
{}

void 
COPKEnergyTicker::initialise()
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      sum[iDim][jDim] = 0.0;
}

void 
COPKEnergyTicker::ticker()
{
  ++count;

  matrix localE;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      localE[iDim][jDim] = 0.0;

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	localE[iDim][jDim] += part.getVelocity()[iDim] * part.getVelocity()[jDim]
	  * Sim->Dynamics.getSpecies(part).getMass();

  //Try and stop round off error this way
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	sum[iDim][jDim] += localE[iDim][jDim];
}

void
COPKEnergyTicker::output(xmlw::XmlStream& XML)
{
  Iflt sumComp(0);
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    sumComp += sum[iDim][iDim];

  XML << xmlw::tag("KEnergyTicker")
      << xmlw::attr("T")
      << sumComp / (((Iflt) count) * ((Iflt)NDIM) * Sim->lN * Sim->Dynamics.units().unitEnergy());
  

  XML << xmlw::tag("KineticTensor");
 
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
      XML << xmlw::tag(name.c_str());

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	  XML << xmlw::attr(name.c_str())
	      << sum[iDim][jDim] / (((Iflt) count) * Sim->lN 				    
				    * Sim->Dynamics.units().unitEnergy());;
	}
      
      XML << xmlw::endtag(name.c_str());
    }
  
  XML << xmlw::endtag("KineticTensor");

  XML << xmlw::endtag("KEnergyTicker");
}

void
COPKEnergyTicker::periodicOutput()
{
  Iflt sumComp(0);
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    sumComp += sum[iDim][iDim];

  sumComp /= ((Iflt) count) * ((Iflt)NDIM) * Sim->lN 
    * Sim->Dynamics.units().unitEnergy();

  I_Pcout() << "<T>_t " <<  sumComp 
	    << ", ";
}

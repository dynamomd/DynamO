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

#include "kenergyticker.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

OPKEnergyTicker::OPKEnergyTicker(const DYNAMO::SimData* tmp, 
		       const XMLNode& XML):
  OPTicker(tmp,"KEnergyTicker"),
  count(0)
{ operator<<(XML); }

void 
OPKEnergyTicker::operator<<(const XMLNode& XML)
{}

void 
OPKEnergyTicker::initialise()
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      sum[iDim][jDim] = 0.0;
}

void 
OPKEnergyTicker::ticker()
{
  ++count;

  matrix localE;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      localE[iDim][jDim] = 0.0;

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	localE[iDim][jDim] += part.getVelocity()[iDim] * part.getVelocity()[jDim]
	  * Sim->dynamics.getSpecies(part).getMass();

  //Try and stop round off error this way
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	sum[iDim][jDim] += localE[iDim][jDim];
}

void
OPKEnergyTicker::output(xml::XmlStream& XML)
{
  Iflt sumComp(0);
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    sumComp += sum[iDim][iDim];

  XML << xml::tag("KEnergyTicker")
      << xml::attr("T")
      << sumComp / (((Iflt) count) * ((Iflt)NDIM) * Sim->N * Sim->dynamics.units().unitEnergy());
  

  XML << xml::tag("KineticTensor");
 
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
      XML << xml::tag(name.c_str());

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	  XML << xml::attr(name.c_str())
	      << sum[iDim][jDim] / (((Iflt) count) * Sim->N 				    
				    * Sim->dynamics.units().unitEnergy());;
	}
      
      XML << xml::endtag(name.c_str());
    }
  
  XML << xml::endtag("KineticTensor");

  XML << xml::endtag("KEnergyTicker");
}

void
OPKEnergyTicker::periodicOutput()
{
  Iflt sumComp(0);
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    sumComp += sum[iDim][iDim];

  sumComp /= ((Iflt) count) * ((Iflt)NDIM) * Sim->N 
    * Sim->dynamics.units().unitEnergy();

  I_Pcout() << "<T>_t " <<  sumComp 
	    << ", ";
}

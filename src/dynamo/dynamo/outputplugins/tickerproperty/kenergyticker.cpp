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

#include <dynamo/outputplugins/tickerproperty/kenergyticker.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPKEnergyTicker::OPKEnergyTicker(const dynamo::SimData* tmp, 
				   const magnet::xml::Node& XML):
    OPTicker(tmp,"KEnergyTicker"),
    count(0)
  { operator<<(XML); }

  void 
  OPKEnergyTicker::operator<<(const magnet::xml::Node& XML)
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
	    * Sim->species[part].getMass(part.getID());

    //Try and stop round off error this way
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	sum[iDim][jDim] += localE[iDim][jDim];
  }

  void
  OPKEnergyTicker::output(magnet::xml::XmlStream& XML)
  {
    double sumComp(0);
  
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      sumComp += sum[iDim][iDim];

    XML << magnet::xml::tag("KEnergyTicker")
	<< magnet::xml::attr("T")
	<< sumComp / (((double) count) * ((double)NDIM) * Sim->N * Sim->units.unitEnergy());
  

    XML << magnet::xml::tag("KineticTensor");
 
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
	XML << magnet::xml::tag(name.c_str());

	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  {
	    std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	    XML << magnet::xml::attr(name.c_str())
		<< sum[iDim][jDim] / (((double) count) * Sim->N 				    
				      * Sim->units.unitEnergy());;
	  }
      
	XML << magnet::xml::endtag(name.c_str());
      }
  
    XML << magnet::xml::endtag("KineticTensor");

    XML << magnet::xml::endtag("KEnergyTicker");
  }

  void
  OPKEnergyTicker::periodicOutput()
  {
    double sumComp(0);
  
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      sumComp += sum[iDim][iDim];

    sumComp /= ((double) count) * ((double)NDIM) * Sim->N 
      * Sim->units.unitEnergy();

    I_Pcout() << "<T>_t " <<  sumComp 
	      << ", ";
  }
}

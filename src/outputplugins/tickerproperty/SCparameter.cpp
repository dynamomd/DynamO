/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "SCparameter.hpp"
#include "../../dynamics/globals/neighbourList.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../dynamics/BC/BC.hpp"
#include <boost/math/special_functions/pow.hpp>  
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <fstream>
#include <cmath>
#include <limits>

OPSCParameter::OPSCParameter(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"SCParameter"),
  maxWaveNumber(0),
  count(0),
  runningsum(0)
{
  operator<<(XML);
}

void 
OPSCParameter::initialise() 
{
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    if (Sim->aspectRatio[iDim] != 1.0) 
      M_throw() << "Cannot use this parameter in a non-cubic box";
  
  maxWaveNumber = lrint(std::pow(Sim->N, 1.0/3.0));

  if (boost::math::pow<3>(maxWaveNumber) != Sim->N)
    M_throw() << "Failed, N does not have an integer cube root!";

  I_cout() << "Max wavelength is "
	   << 1.0 / (maxWaveNumber * Sim->dynamics.units().unitLength());

  maxWaveNumber *= 2;

  runningsum.resize(maxWaveNumber + 1, 0);

  ticker();
}

void 
OPSCParameter::ticker()
{
  ++count;

  for (size_t k(0); k <= maxWaveNumber; ++k)
    {
      std::complex<double> sum(0, 0);

      BOOST_FOREACH(const Particle& part, Sim->particleList)
	{
	  double psum(0);
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    psum += part.getPosition()[iDim];
	  
	  psum *= 2.0 * M_PI * k;
	  sum += std::complex<double>(std::cos(psum), std::sin(psum));
	}
      
      runningsum[k] += std::abs(sum);
    }
}

void 
OPSCParameter::output(xml::XmlStream& XML)
{
  XML << xml::tag("SCParameter")
      << xml::attr("SCWaveNumber") 
      << lrint(std::pow(Sim->N, 1.0/3.0))
      << xml::attr("SCWaveNumberVal") 
      << runningsum[lrint(std::pow(Sim->N, 1.0/3.0))] 
    / (static_cast<double>(count) * Sim->N) 
      << xml::chardata();
  
  for (size_t k(0); k <= maxWaveNumber; ++k)
    {
      XML << k * Sim->dynamics.units().unitLength() << " "
	  << runningsum[k] / (static_cast<double>(count) * Sim->N) 
	  << "\n";
    }

  XML << xml::endtag("SCParameter");
}

void 
OPSCParameter::operator<<(const XMLNode& XML)
{
}

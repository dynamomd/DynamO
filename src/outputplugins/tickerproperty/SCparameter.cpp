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

#include "SCparameter.hpp"
#include <fstream>
#include <cmath>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include <limits>
#include "../../dynamics/globals/neighbourList.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../dynamics/BC/BC.hpp"
#include <boost/math/special_functions/pow.hpp>  

COPSCParameter::COPSCParameter(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPTicker(tmp,"SCParameter"),
  maxWaveNumber(0),
  count(0),
  runningsum(0)
{
  operator<<(XML);
}

void 
COPSCParameter::initialise() 
{
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    if (Sim->aspectRatio[iDim] != 1.0) 
      D_throw() << "Cannot use this parameter in a non-cubic box";
  
  maxWaveNumber = lrint(std::pow(Sim->lN, 1.0/3.0));

  if (boost::math::pow<3>(maxWaveNumber) != Sim->lN)
    D_throw() << "Failed, N does not have an integer cube root!";

  I_cout() << "Max wavelength is "
	   << 1.0 / (maxWaveNumber * Sim->Dynamics.units().unitLength());

  runningsum.resize(maxWaveNumber + 1, 0);

  ticker();
}

void 
COPSCParameter::ticker()
{
  ++count;

  for (size_t k(0); k <= maxWaveNumber; ++k)
    {
      std::complex<Iflt> sum(0, 0);

      BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
	{
	  Iflt psum(0);
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    psum += part.getPosition()[iDim];
	  
	  psum *= 2.0 * PI * k;
	  sum += std::complex<Iflt>(std::cos(psum), std::sin(psum));
	}
      
      runningsum[k] += std::abs(sum);
    }
}

void 
COPSCParameter::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("SCParameter")
      << xmlw::attr("MaxWaveVal") 
      << runningsum.back() / (static_cast<Iflt>(count) * Sim->lN)
      << xmlw::chardata();
  
  for (size_t k(0); k <= maxWaveNumber; ++k)
    {
      XML << k * Sim->Dynamics.units().unitLength() << " "
	  << runningsum[k] / (static_cast<Iflt>(count) * Sim->lN) 
	  << "\n";
    }

  XML << xmlw::endtag("SCParameter");
}

void 
COPSCParameter::operator<<(const XMLNode& XML)
{
}

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

#include "SHcrystal.hpp"
#include <fstream>
#include <cmath>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include <limits>
#include "../../dynamics/globals/neighbourList.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../dynamics/BC/BC.hpp"
#include <boost/math/special_functions/spherical_harmonic.hpp>

COPSHCrystal::COPSHCrystal(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"SHCrystal"), rg(1), maxl(6),
  nblistID(std::numeric_limits<size_t>::max()),
  count(0)
{}

void 
COPSHCrystal::initialise() 
{ 
  rg = 1.5 * Sim->Dynamics.units().unitLength();

  double smallestlength = HUGE_VAL;
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& pGlob, Sim->Dynamics.getGlobals())
    if (dynamic_cast<const CGNeighbourList*>(pGlob.get_ptr()) != NULL)
      {
	const Iflt l(static_cast<const CGNeighbourList*>(pGlob.get_ptr())
		     ->getMaxSupportedInteractionLength());
	if ((l >= rg) && (l < smallestlength))
	  {
	    //this neighbourlist is better suited
	    smallestlength = l;
	    nblistID = pGlob->getID();
	  }
      }

  if (nblistID == std::numeric_limits<size_t>::max())
    D_throw() << "There is not a suitable neighbourlist for the cut-off radius selected."
      "\nR_g = " << rg / Sim->Dynamics.units().unitLength();

  globalcoeff.resize(maxl);
  for (size_t l=0; l < maxl; ++l)
    globalcoeff[l].resize(2*l+1,std::complex<Iflt>(0,0));

  ticker();
}
void 
COPSHCrystal::ticker()
{
  ++count;
  sphericalsum ssum(Sim, rg, maxl);
  
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      static_cast<const CGNeighbourList*>
	(Sim->Dynamics.getGlobals()[nblistID].get_ptr())
	->getParticleNeighbourhood
	(part, CGNeighbourList::getNBDelegate
	 (&sphericalsum::operator(), &ssum));
      
      for (size_t l(0); l < maxl; ++l)
	for (int m(-l); m <= static_cast<int>(l); ++m)
	  globalcoeff[l][m+l] += ssum.coeffsum[l][m+l] 
	    / std::complex<Iflt>(ssum.count,0);
      
      ssum.clear();
    }

}

void 
COPSHCrystal::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("SHCrystal");
  
  for (size_t l(0); l < maxl; ++l)
    {
      XML << xmlw::tag("Q")
	  << xmlw::attr("l") << l;

      Iflt sum(0);
      for (int m(-l); m <= static_cast<int>(l); ++m)
	sum += std::norm(globalcoeff[l][m+l]);
      
      XML << xmlw::attr("val") 
	  << sum * 4.0 * PI / (2.0 * l + 1.0)
	  << xmlw::endtag("Q");
    }

  XML << xmlw::endtag("SHCrystal");
}

COPSHCrystal::sphericalsum::sphericalsum
(const DYNAMO::SimData * const nSim, const Iflt& nrg, const size_t& nl):
  Sim(nSim), rg(nrg), maxl(nl), count(0)
{
  coeffsum.resize(maxl);
  for (size_t l=0; l < maxl; ++l)
    coeffsum[l].resize(2*l+1,std::complex<Iflt>(0,0));
}

void 
COPSHCrystal::sphericalsum::operator()
  (const CParticle& part, const size_t& ID) const
{
  Vector rij = part.getPosition() - Sim->vParticleList[ID].getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  
  Iflt norm = rij.nrm();
  if (norm <= rg)
    {
      ++count;
      rij /= norm;
      Iflt theta = std::acos(rij[0]);
      Iflt phi = std::asin(rij[1] / std::sin(theta));

      for (size_t l(0); l < maxl; ++l)
	for (int m(-l); m <= static_cast<int>(l); ++m)
	  coeffsum[l][m+l] += boost::math::spherical_harmonic(l,m,theta,phi);
    }
}

void
COPSHCrystal::sphericalsum::clear()
{
  for (size_t l(0); l < maxl; ++l)
    for (int m(-l); m <= static_cast<int>(l); ++m)
      coeffsum[l][m+l] = std::complex<Iflt>(0,0);
}

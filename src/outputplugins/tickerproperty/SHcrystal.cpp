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
#include "../../extcode/wignerThreeJ.hpp"

COPSHCrystal::COPSHCrystal(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPTicker(tmp,"SHCrystal"), rg(1.2), maxl(7),
  nblistID(std::numeric_limits<size_t>::max()),
  count(0)
{
  operator<<(XML);
}

void 
COPSHCrystal::initialise() 
{ 
  

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
	  globalcoeff[l][m+l] += ssum.coeffsum[l][m+l];
      
      count += ssum.count;

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
      
      Iflt Qsum(0);
      for (int m(-l); m <= static_cast<int>(l); ++m)
	Qsum += std::norm(globalcoeff[l][m+l] / std::complex<Iflt>(count, 0));
      
      XML << xmlw::attr("val")
	  << std::sqrt(Qsum * 4.0 * PI / (2.0 * l + 1.0))
	  << xmlw::endtag("Q");
      
      XML << xmlw::tag("W")
	  << xmlw::attr("l") << l;

      std::complex<Iflt> Wsum(0,0);
      for (int m1(-l); m1 <= static_cast<int>(l); ++m1)
	for (int m2(-l); m2 <= static_cast<int>(l); ++m2)
	  if (std::abs(m1 + m2) <= static_cast<int>(l))
	    Wsum += DYNAMO::threej(l,l,l,m1,m2,-(m1+m2))
	      * globalcoeff[l][m1] * globalcoeff[l][m2]
	      * globalcoeff[l][-(m1+m2)] 
	      / std::complex<Iflt>(count * count * count, 0);
      
      XML << xmlw::attr("val")
	  << Wsum * std::pow(Qsum, -1.5)
	  << xmlw::endtag("W");
    }

  XML << xmlw::endtag("SHCrystal");
}

void 
COPSHCrystal::operator<<(const XMLNode& XML)
{
  rg *= Sim->Dynamics.units().unitLength();

  try
    {
      if (XML.isAttributeSet("CutOffR"))
	rg = boost::lexical_cast<Iflt>(XML.getAttribute("CutOffR"))
	  * Sim->Dynamics.units().unitLength();

      if (XML.isAttributeSet("MaxL"))
	maxl = boost::lexical_cast<size_t>(XML.getAttribute("MaxL"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPSHCrystal";
    }

  I_cout() << "Cut off radius of " 
	   << rg / Sim->Dynamics.units().unitLength();
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
      Iflt sintheta = std::sin(theta);
      Iflt phi = rij[1] / sintheta;
      
      if (fabs(phi) > 1.0)
	phi = (phi > 0) ? 0.5 * PI : 1.5 * PI;
      else
	phi = std::asin(phi);

      if (std::sin(theta) == 0) phi = 0;
      
      if (phi < 0) phi += 2.0*PI;

      for (size_t l(0); l < maxl; ++l)
	for (int m(-l); m <= static_cast<int>(l); ++m)
	  coeffsum[l][m+l] += boost::math::spherical_harmonic(l,m,theta,phi);
    }
}

void
COPSHCrystal::sphericalsum::clear()
{
  count = 0;

  for (size_t l(0); l < maxl; ++l)
    for (int m(-l); m <= static_cast<int>(l); ++m)
      coeffsum[l][m+l] = std::complex<Iflt>(0,0);
}


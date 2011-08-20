/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/outputplugins/tickerproperty/SHcrystal.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/foreach.hpp>
#include <magnet/math/wigner3J.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>
#include <cmath>
#include <limits>

namespace dynamo {
  OPSHCrystal::OPSHCrystal(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"SHCrystal"), rg(1.2), maxl(7),
    nblistID(std::numeric_limits<size_t>::max()),
    count(0)
  {
    operator<<(XML);
  }

  void 
  OPSHCrystal::operator<<(const magnet::xml::Node& XML)
  {
    try
      {
	if (XML.hasAttribute("CutOffR"))
	  rg = XML.getAttribute("CutOffR").as<double>();

	if (XML.hasAttribute("MaxL"))
	  maxl = XML.getAttribute("MaxL").as<size_t>();
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPSHCrystal";
      }

    rg *= Sim->dynamics.units().unitLength();


    dout << "Cut off radius of " 
	 << rg / Sim->dynamics.units().unitLength() << std::endl;
  }


  void 
  OPSHCrystal::initialise() 
  { 
    double smallestlength = HUGE_VAL;
    BOOST_FOREACH(const magnet::ClonePtr<Global>& pGlob, Sim->dynamics.getGlobals())
      if (dynamic_cast<const GNeighbourList*>(pGlob.get_ptr()) != NULL)
	{
	  const double l(static_cast<const GNeighbourList*>(pGlob.get_ptr())
			 ->getMaxSupportedInteractionLength());
	  if ((l >= rg) && (l < smallestlength))
	    {
	      //this neighbourlist is better suited
	      smallestlength = l;
	      nblistID = pGlob->getID();
	    }
	}

    if (nblistID == std::numeric_limits<size_t>::max())
      M_throw() << "There is not a suitable neighbourlist for the cut-off radius selected."
	"\nR_g = " << rg / Sim->dynamics.units().unitLength();

    globalcoeff.resize(maxl);
    for (size_t l=0; l < maxl; ++l)
      globalcoeff[l].resize(2*l+1,std::complex<double>(0,0));

    ticker();
  }

  void 
  OPSHCrystal::ticker()
  {
    sphericalsum ssum(Sim, rg, maxl);
  
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	static_cast<const GNeighbourList*>
	  (Sim->dynamics.getGlobals()[nblistID].get_ptr())
	  ->getParticleNeighbourhood
	  (part, magnet::function::MakeDelegate(&ssum, &sphericalsum::operator()));
      
	for (size_t l(0); l < maxl; ++l)
	  for (int m(-l); m <= static_cast<int>(l); ++m)
	    globalcoeff[l][m+l] += ssum.coeffsum[l][m+l];
      
	count += ssum.count;

	ssum.clear();
      }
  }

  void 
  OPSHCrystal::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("SHCrystal");
  
    for (int l(0); l < static_cast<int>(maxl); ++l)
      {
	XML << magnet::xml::tag("Q")
	    << magnet::xml::attr("l") << l;
      
	double Qsum(0);
	for (int m(-l); m <= l; ++m)
	  Qsum += std::norm(globalcoeff[l][m+l] / std::complex<double>(count, 0));
      
	XML << magnet::xml::attr("val")
	    << std::sqrt(Qsum * 4.0 * M_PI / (2.0 * l + 1.0))
	    << magnet::xml::endtag("Q");
      
	XML << magnet::xml::tag("W")
	    << magnet::xml::attr("l") << l;

	std::complex<double> Wsum(0, 0);
	for (int m1(-l); m1 <= l; ++m1)
	  for (int m2(-l); m2 <= l; ++m2)
	    {
	      int m3 = -(m1+m2);
	      if (std::abs(m3) <= l)
		Wsum += std::complex<double>(magnet::math::wignerThreej(l,l,l,m1,m2,m3) 
					     * std::pow(count,-3.0), 0)
		  * globalcoeff[l][m1+l]
		  * globalcoeff[l][m2+l]
		  * globalcoeff[l][l-(m1+m2)]
		  ;
	    }
      
	XML << magnet::xml::attr("val")
	    << Wsum * std::pow(Qsum, -1.5)
	    << magnet::xml::endtag("W");
      }

    XML << magnet::xml::endtag("SHCrystal");
  }

  OPSHCrystal::sphericalsum::sphericalsum
  (const dynamo::SimData * const nSim, const double& nrg, const size_t& nl):
    Sim(nSim), rg(nrg), maxl(nl), count(0)
  {
    coeffsum.resize(maxl);
    for (size_t l=0; l < maxl; ++l)
      coeffsum[l].resize(2*l+1,std::complex<double>(0,0));
  }

  void 
  OPSHCrystal::sphericalsum::operator()
    (const Particle& part, const size_t& ID) const
  {
    if (part.getID() == ID) return;
    Vector rij = part.getPosition() - Sim->particleList[ID].getPosition();
    Sim->dynamics.BCs().applyBC(rij);
  
    double norm = rij.nrm();
    if (norm <= rg)
      {
	++count;
	rij /= norm;
	double theta = std::acos(rij[0]);
	double sintheta = std::sin(theta);
	double phi = rij[1] / sintheta;
      
	if (fabs(phi) > 1.0)
	  phi = (phi > 0) ? 0.5 * M_PI : 1.5 * M_PI;
	else
	  phi = std::asin(phi);

	if (std::sin(theta) == 0) phi = 0;
      
	if (phi < 0) phi += 2.0 * M_PI;

	for (size_t l(0); l < maxl; ++l)
	  for (int m(-l); m <= static_cast<int>(l); ++m)
	    coeffsum[l][m+l] += boost::math::spherical_harmonic(l,m,theta,phi);
      }
  }

  void
  OPSHCrystal::sphericalsum::clear()
  {
    count = 0;

    for (size_t l(0); l < maxl; ++l)
      for (int m(-l); m <= static_cast<int>(l); ++m)
	coeffsum[l][m+l] = std::complex<double>(0,0);
  }
}

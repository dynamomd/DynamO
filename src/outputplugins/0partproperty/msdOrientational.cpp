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

#include "msdOrientational.hpp"
#include <boost/foreach.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include <vector>


OPMSDOrientational::OPMSDOrientational(const DYNAMO::SimData* tmp, const XMLNode&):
  OutputPlugin(tmp,"MSDOrientational")
{}

OPMSDOrientational::~OPMSDOrientational()
{}

void
OPMSDOrientational::initialise()
{
  if (dynamic_cast<const LNOrientation*>(&Sim->dynamics.getLiouvillean()) == NULL)
  {
    D_throw() << "Plugin requires species to define an orientation";
  }

  initialConfiguration.clear();
  initialConfiguration.resize(Sim->N);

  const std::vector<LNOrientation::rotData>& rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getCompleteRotData());

  for (size_t ID = 0; ID < Sim->N; ++ID)
  {
    initialConfiguration[ID] = RUpair(Sim->particleList[ID].getPosition(), rdat[ID].orientation);
  }
}

void
OPMSDOrientational::output(xml::XmlStream &XML)
{
  msdCalcReturn MSDOrientational(calculate());

  XML << xml::tag("MSDOrientational")

      << xml::tag("Perpendicular")
      << xml::attr("val") << MSDOrientational.perpendicular
      << xml::attr("diffusionCoeff")
      << MSDOrientational.perpendicular * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xml::endtag("Perpendicular")

      << xml::tag("Parallel")
      << xml::attr("val") << MSDOrientational.parallel
      << xml::attr("diffusionCoeff")
      << MSDOrientational.parallel * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xml::endtag("Parallel")

      << xml::tag("Rotational")
      << xml::attr("method") << "LegendrePolynomial1"
      << xml::attr("val") << MSDOrientational.rotational_legendre1
      << xml::attr("diffusionCoeff")
      << MSDOrientational.rotational_legendre1 * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xml::endtag("Rotational")

      << xml::tag("Rotational")
      << xml::attr("method") << "LegendrePolynomial2"
      << xml::attr("val") << MSDOrientational.rotational_legendre2
      << xml::attr("diffusionCoeff")
      << MSDOrientational.rotational_legendre2 * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xml::endtag("Rotational")

      << xml::endtag("MSDOrientational");
}

OPMSDOrientational::msdCalcReturn
OPMSDOrientational::calculate() const
{
  msdCalcReturn MSR;

  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();

  Iflt 	acc_perp(0.0), acc_parallel(0.0), longitudinal_projection(0.0),
	acc_rotational_legendre1(0.0), acc_rotational_legendre2(0.0),
	cos_theta(0.0);

  Vector displacement_term(0,0,0);

  const std::vector<LNOrientation::rotData>& latest_rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getCompleteRotData());

  BOOST_FOREACH(const Particle& part, Sim->particleList)
  {
    displacement_term = part.getPosition() - initialConfiguration[part.getID()].first;
    longitudinal_projection = (displacement_term | initialConfiguration[part.getID()].second);
    cos_theta = (initialConfiguration[part.getID()].second | latest_rdat[part.getID()].orientation);

    acc_perp += (displacement_term - (longitudinal_projection * initialConfiguration[part.getID()].second)).nrm2();
    acc_parallel += std::pow(longitudinal_projection, 2);
    acc_rotational_legendre1 += boost::math::legendre_p(1, cos_theta);
    acc_rotational_legendre2 += boost::math::legendre_p(2, cos_theta);
  }

  // In the N-dimensional case, the parallel component is 1-dimensonal and the perpendicular (N-1)
  acc_perp /= (initialConfiguration.size() * 2.0 * (NDIM - 1) * Sim->dynamics.units().unitArea());
  acc_parallel /= (initialConfiguration.size() * 2.0 * Sim->dynamics.units().unitArea());

  // Rotational forms by Magda, Davis and Tirrell
  // <P1(cos(theta))> = exp[-2 D t]
  // <P2(cos(theta))> = exp[-6 D t]

  // WARNING! - only valid for sufficiently high density
  // Use the MSDOrientationalCorrelator to check for an exponential fit

  acc_rotational_legendre1 = log(acc_rotational_legendre1 / initialConfiguration.size()) / (-2.0);
  acc_rotational_legendre2 = log(acc_rotational_legendre2 / initialConfiguration.size()) / (-6.0);

  MSR.parallel = acc_parallel;
  MSR.perpendicular = acc_perp;
  MSR.rotational_legendre1 = acc_rotational_legendre1;
  MSR.rotational_legendre2 = acc_rotational_legendre2;

  return MSR;
}

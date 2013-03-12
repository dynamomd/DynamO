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

#include <dynamo/outputplugins/0partproperty/msdOrientational.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
  OPMSDOrientational::OPMSDOrientational(const dynamo::Simulation* tmp, 
					 const magnet::xml::Node&):
    OutputPlugin(tmp,"MSDOrientational")
  {}

  OPMSDOrientational::~OPMSDOrientational()
  {}

  void
  OPMSDOrientational::initialise()
  {
    initialConfiguration.clear();
    initialConfiguration.resize(Sim->N);

    const std::vector<Dynamics::rotData>& rdat(Sim->dynamics->getCompleteRotData());

    for (size_t ID = 0; ID < Sim->N; ++ID)
      {
	initialConfiguration[ID] = RUpair(Sim->particles[ID].getPosition(), rdat[ID].orientation * Quaternion::initialDirector());
      }
  }

  void
  OPMSDOrientational::output(magnet::xml::XmlStream &XML)
  {
    msdCalcReturn MSDOrientational(calculate());

    XML << magnet::xml::tag("MSDOrientational")

	<< magnet::xml::tag("Perpendicular")
	<< magnet::xml::attr("val") << MSDOrientational.perpendicular
	<< magnet::xml::attr("diffusionCoeff")
	<< MSDOrientational.perpendicular * Sim->units.unitTime() / Sim->systemTime
	<< magnet::xml::endtag("Perpendicular")

	<< magnet::xml::tag("Parallel")
	<< magnet::xml::attr("val") << MSDOrientational.parallel
	<< magnet::xml::attr("diffusionCoeff")
	<< MSDOrientational.parallel * Sim->units.unitTime() / Sim->systemTime
	<< magnet::xml::endtag("Parallel")

	<< magnet::xml::tag("Rotational")
	<< magnet::xml::attr("method") << "LegendrePolynomial1"
	<< magnet::xml::attr("val") << MSDOrientational.rotational_legendre1
	<< magnet::xml::attr("diffusionCoeff")
	<< MSDOrientational.rotational_legendre1 * Sim->units.unitTime() / Sim->systemTime
	<< magnet::xml::endtag("Rotational")

	<< magnet::xml::tag("Rotational")
	<< magnet::xml::attr("method") << "LegendrePolynomial2"
	<< magnet::xml::attr("val") << MSDOrientational.rotational_legendre2
	<< magnet::xml::attr("diffusionCoeff")
	<< MSDOrientational.rotational_legendre2 * Sim->units.unitTime() / Sim->systemTime
	<< magnet::xml::endtag("Rotational")

	<< magnet::xml::endtag("MSDOrientational");
  }

  OPMSDOrientational::msdCalcReturn
  OPMSDOrientational::calculate() const
  {
    msdCalcReturn MSR;

    //Required to get the correct results
    Sim->dynamics->updateAllParticles();

    double 	acc_perp(0.0), acc_parallel(0.0), longitudinal_projection(0.0),
      acc_rotational_legendre1(0.0), acc_rotational_legendre2(0.0),
      cos_theta(0.0);

    Vector displacement_term(0,0,0);

    const std::vector<Dynamics::rotData>& latest_rdat(Sim->dynamics->getCompleteRotData());

    BOOST_FOREACH(const Particle& part, Sim->particles)
      {
	displacement_term = part.getPosition() - initialConfiguration[part.getID()].first;
	longitudinal_projection = (displacement_term | initialConfiguration[part.getID()].second);
	cos_theta = (initialConfiguration[part.getID()].second | (latest_rdat[part.getID()].orientation * Quaternion::initialDirector()));

	acc_perp += (displacement_term - (longitudinal_projection * initialConfiguration[part.getID()].second)).nrm2();
	acc_parallel += std::pow(longitudinal_projection, 2);
	double clamp_cos_theta = std::max(std::min(1.0, cos_theta), -1.0);
	acc_rotational_legendre1 += boost::math::legendre_p(1, clamp_cos_theta);
	acc_rotational_legendre2 += boost::math::legendre_p(2, clamp_cos_theta);
      }

    // In the N-dimensional case, the parallel component is 1-dimensonal and the perpendicular (N-1)
    acc_perp /= (initialConfiguration.size() * 2.0 * (NDIM - 1) * Sim->units.unitArea());
    acc_parallel /= (initialConfiguration.size() * 2.0 * Sim->units.unitArea());

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
}

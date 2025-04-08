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

#include <boost/math/special_functions/legendre.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/include.hpp>
#include <dynamo/outputplugins/tickerproperty/msdOrientationalCorrelator.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/sysTicker.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPMSDOrientationalCorrelator::OPMSDOrientationalCorrelator(
    const dynamo::Simulation *tmp, const magnet::xml::Node &XML)
    : OPTicker(tmp, "MSDOrientationalCorrelator"), length(50),
      currCorrLength(0), ticksTaken(0), notReady(true) {
  OPMSDOrientationalCorrelator::operator<<(XML);
}

void OPMSDOrientationalCorrelator::operator<<(const magnet::xml::Node &XML) {
  if (XML.hasAttribute("Length"))
    length = XML.getAttribute("Length").as<size_t>();
}

void OPMSDOrientationalCorrelator::initialise() {
  dout << "The length of the MSD orientational correlator is " << length
       << std::endl;

  historicalData.resize(Sim->N(), boost::circular_buffer<RUpair>(length));

  stepped_data_parallel.resize(length, double(0.0));
  stepped_data_perpendicular.resize(length, double(0.0));
  stepped_data_rotational_legendre1.resize(length, double(0.0));
  stepped_data_rotational_legendre2.resize(length, double(0.0));

  // The Legendre polynomials are equal to 1 at t = 0
  stepped_data_rotational_legendre1[0] = 1.0;
  stepped_data_rotational_legendre2[0] = 1.0;

  currCorrLength = 1.0;

  const std::vector<Dynamics::rotData> &initial_rdat(
      Sim->dynamics->getCompleteRotData());

  for (const Particle &part : Sim->particles) {
    historicalData[part.getID()].push_front(
        RUpair(part.getPosition(), initial_rdat[part.getID()].orientation *
                                       Quaternion::initialDirector()));
  }
}

void OPMSDOrientationalCorrelator::ticker() {
  const std::vector<Dynamics::rotData> &current_rdat(
      Sim->dynamics->getCompleteRotData());
  for (const Particle &part : Sim->particles) {
    historicalData[part.getID()].push_front(
        RUpair(part.getPosition(), current_rdat[part.getID()].orientation *
                                       Quaternion::initialDirector()));
  }

  if (notReady) {
    if (++currCorrLength != length) {
      return;
    }

    notReady = false;
  }

  accPass();
}

void OPMSDOrientationalCorrelator::accPass() {
  ++ticksTaken;

  double longitudinal_projection(0.0), cos_theta(0.0);
  Vector displacement_term{0, 0, 0};

  for (const Particle &part : Sim->particles) {
    for (size_t step(0); step < length; ++step) {
      displacement_term = historicalData[part.getID()][step].first -
                          historicalData[part.getID()][0].first;
      longitudinal_projection =
          (displacement_term | historicalData[part.getID()][0].second);
      cos_theta = (historicalData[part.getID()][step].second |
                   historicalData[part.getID()][0].second);

      stepped_data_parallel[step] += std::pow(longitudinal_projection, 2);
      stepped_data_perpendicular[step] +=
          (displacement_term -
           (longitudinal_projection * historicalData[part.getID()][0].second))
              .nrm2();

      double clamp_cos_theta = std::max(std::min(1.0, cos_theta), -1.0);
      stepped_data_rotational_legendre1[step] +=
          boost::math::legendre_p(1, clamp_cos_theta);
      stepped_data_rotational_legendre2[step] +=
          boost::math::legendre_p(2, clamp_cos_theta);
    }
  }
}

void OPMSDOrientationalCorrelator::output(magnet::xml::XmlStream &XML) {
  // Begin XML output
  XML << magnet::xml::tag("MSDOrientationalCorrelator");

  double dt = dynamic_cast<const SysTicker &>(*Sim->systems["SystemTicker"])
                  .getPeriod() /
              Sim->units.unitTime();

  XML << magnet::xml::tag("Component") << magnet::xml::attr("Type")
      << "Parallel" << magnet::xml::chardata();

  for (size_t step(0); step < length; ++step) {
    XML << dt * step << "\t"
        << stepped_data_parallel[step] /
               (static_cast<double>(ticksTaken) *
                static_cast<double>(Sim->N()) * Sim->units.unitArea())
        << "\n";
  }

  XML << magnet::xml::endtag("Component");

  XML << magnet::xml::tag("Component") << magnet::xml::attr("Type")
      << "Perpendicular" << magnet::xml::chardata();

  for (size_t step(0); step < length; ++step) {
    XML << dt * step << "\t"
        << stepped_data_perpendicular[step] /
               (static_cast<double>(ticksTaken) *
                static_cast<double>(Sim->N()) * Sim->units.unitArea())
        << "\n";
  }

  XML << magnet::xml::endtag("Component");

  XML << magnet::xml::tag("Component") << magnet::xml::attr("Type")
      << "Rotational";

  XML << magnet::xml::tag("Method") << magnet::xml::attr("Name")
      << "LegendrePolynomial1" << magnet::xml::chardata();

  for (size_t step(0); step < length; ++step) {
    XML << dt * step << "\t"
        << stepped_data_rotational_legendre1[step] /
               (static_cast<double>(ticksTaken) * static_cast<double>(Sim->N()))
        << "\n";
  }

  XML << magnet::xml::endtag("Method");

  XML << magnet::xml::tag("Method") << magnet::xml::attr("Name")
      << "LegendrePolynomial2" << magnet::xml::chardata();

  for (size_t step(0); step < length; ++step) {
    XML << dt * step << "\t"
        << stepped_data_rotational_legendre2[step] /
               (static_cast<double>(ticksTaken) * static_cast<double>(Sim->N()))
        << "\n";
  }

  XML << magnet::xml::endtag("Method");

  XML << magnet::xml::endtag("Component");

  XML << magnet::xml::endtag("MSDOrientationalCorrelator");
}
} // namespace dynamo

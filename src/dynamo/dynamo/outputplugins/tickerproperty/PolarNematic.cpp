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

#include <cmath>
#include <complex>
#include <dynamo/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/tickerproperty/PolarNematic.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>
#include <fstream>
#include <limits>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPPolarNematic::OPPolarNematic(const dynamo::Simulation *tmp,
                               const magnet::xml::Node &XML)
    : OPTicker(tmp, "PolarNematic") {
  operator<<(XML);
}

void OPPolarNematic::operator<<(const magnet::xml::Node &XML) {}

void OPPolarNematic::initialise() { ticker(); }

void OPPolarNematic::ticker() {
  std::complex<double> polar(0, 0);
  std::complex<double> nematic(0, 0);
  const auto &data = Sim->dynamics->getCompleteRotData();
  for (const auto &entry : data) {
    const Vector vec = entry.orientation * Quaternion::initialDirector();
    const double arg = std::arg(std::complex<double>(vec[0], vec[1]));
    polar += std::exp(std::complex<double>(0, arg));
    nematic += std::exp(std::complex<double>(0, 2 * arg));
  }

  const size_t count = data.size();
  polar /= count;
  nematic /= count;

  _history.push_back(
      std::pair<double, double>(std::abs(polar), std::abs(nematic)));
}

void OPPolarNematic::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("PolarNematic");

  double polar(0), nematic(0);
  for (const auto &val : _history) {
    polar += val.first;
    nematic += val.second;
  }

  XML << magnet::xml::attr("PolarAvg") << polar / _history.size()
      << magnet::xml::attr("NematicAvg") << nematic / _history.size()
      << magnet::xml::chardata();

  for (const auto &val : _history)
    XML << "\n" << val.first << " " << val.second;

  XML << magnet::xml::endtag("PolarNematic");
}
} // namespace dynamo

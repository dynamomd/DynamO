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

#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

enum EnergyMode {
  MIDPOINT,
  LEFT,
  RIGHT,
  VOLUME
};

namespace dynamo {
  void 
  PotentialLennardJones::outputXML(magnet::xml::XmlStream& XML) const {
    XML << magnet::xml::attr("Type") << "LennardJones"
	<< magnet::xml::attr("CutOff") << _cutoff
	<< magnet::xml::attr("DeltaE") << _deltaE
      ;
  }

  double 
  PotentialLennardJones::U_uncut(double r) const {
    return 4 * _epsilon * (std::pow(_sigma / r, 12.0) - std::pow(_sigma / r, 6.0));
  }

  double 
  PotentialLennardJones::U(double r) const {
    return U(r) - U(_cutoff);
  }

  double 
  PotentialLennardJones::minimum() const {
    return _sigma * std::pow(2.0, 1.0/6.0);
  }

  void 
  PotentialLennardJones::operator<<(const magnet::xml::Node& XML) {
    _cutoff = XML.getAttribute("CutOff").as<double>();
    _deltaE = XML.getAttribute("DeltaE").as<double>();
    _sigma = XML.getAttribute("Sigma").as<double>();
    _epsilon = XML.getAttribute("Epsilon").as<double>();
    
    std::string mode_string = XML.getAttribute("Mode").as<std::string>();
    if (!mode_string.compare("Midpoint"))    _mode = MIDPOINT;
    else if (!mode_string.compare("Left"))   _mode = LEFT;
    else if (!mode_string.compare("Right"))  _mode = RIGHT;
    else if (!mode_string.compare("Volume")) _mode = VOLUME;
    else
      M_throw() << "Unknown LennardJones Mode (" << mode_string << ") at " << XML.getPath();
  }

  void 
  PotentialLennardJones::calculateToStep(size_t id) const {
    M_throw() << "Not implemented";
  }
}

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

#include <dynamo/interactions/potentials/potential.hpp>
#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  shared_ptr<Potential>
  Potential::getClass(const magnet::xml::Node& XML)
  {
    if (!XML.getAttribute("Type").getValue().compare("Stepped"))
      return shared_ptr<Potential>(new PotentialStepped(XML));
    else if (!XML.getAttribute("Type").getValue().compare("LennardJones"))
      return shared_ptr<Potential>(new PotentialLennardJones(XML));
    else 
      M_throw() << XML.getAttribute("Type").getValue()
		<< ", Unknown type of Potential encountered";
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const Potential& g)
  {
    XML << magnet::xml::tag("Potential");
    g.outputXML(XML);
    XML << magnet::xml::endtag("Potential");
    return XML;
  }

  PotentialStepped::PotentialStepped(std::vector<std::pair<double, double> > steps, bool direction):
    _direction(direction)
  {
    std::sort(steps.rbegin(), steps.rend());
    for (size_t i(0); i < steps.size(); ++i)
      {
	_r_cache.push_back(steps[i].first);
	_u_cache.push_back(steps[i].second);
      }
  }

  void 
  PotentialStepped::outputXML(magnet::xml::XmlStream& XML) const {
    XML << magnet::xml::attr("Type") << "Stepped";

    if (_direction)
      XML << magnet::xml::attr("Direction") << "Right";
    else
      XML << magnet::xml::attr("Direction") << "Left";

    for (size_t id(0); id < steps(); ++id)
      XML << magnet::xml::tag("Step")
	  << magnet::xml::attr("R") << _r_cache[id]
	  << magnet::xml::attr("E") << _u_cache[id]
	  << magnet::xml::endtag("Step");
  }

  void 
  PotentialStepped::operator<<(const magnet::xml::Node& XML) {
    if (XML.getAttribute("Direction").as<std::string>() == "Left")
      _direction = false;
    else if (XML.getAttribute("Direction").as<std::string>() == "Right")
      _direction = true;
    else
      M_throw() << "Could not parse Direction, should be either \"Left\" or \"Right\"\nXML path:"
		<< XML.getPath();

    std::vector<std::pair<double, double> > steps;
    for (magnet::xml::Node node = XML.fastGetNode("Step"); node.valid(); ++node)
      steps.push_back(std::pair<double, double>(node.getAttribute("R").as<double>(),
						node.getAttribute("E").as<double>()));
    
    if (steps.empty())
      M_throw() << "You cannot load a stepped potential with no steps.\nXML path: " << XML.getPath();
    
    *this = PotentialStepped(steps, _direction);
  }
}

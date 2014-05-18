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

#include <dynamo/outputplugins/tickerproperty/PolarNematic.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>
#include <cmath>
#include <limits>

namespace dynamo {
  OPPolarNematic::OPPolarNematic(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"PolarNematic")
  {
    operator<<(XML);
  }

  void 
  OPPolarNematic::operator<<(const magnet::xml::Node& XML)
  {}


  void 
  OPPolarNematic::initialise() 
  { ticker(); }

  void 
  OPPolarNematic::ticker()
  {
    Vector polar(0,0,0);
    Vector nematic(0,0,0);
    for (const auto& data : Sim->dynamics->getCompleteRotData())
      {
	Vector vec = data.orientation * Quaternion::initialDirector();
	polar += vec;
	for (size_t i(0); i < 3; ++i)
	  nematic[i] += 2 * vec[i] * vec[i] - 1;
      }

    const size_t count = Sim->dynamics->getCompleteRotData().size();
    polar /= count;
    nematic /= count;

    _history.push_back(std::pair<double, double>(polar.nrm(), nematic.nrm()));
  }

  void 
  OPPolarNematic::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("PolarNematic");

    double polar(0), nematic(0);
    for (const auto& val : _history)
      { polar += val.first; nematic += val.second; }
    
    XML << magnet::xml::attr("PolarAvg") << polar / _history.size()
	<< magnet::xml::attr("NematicAvg") << nematic / _history.size()
	<< magnet::xml::chardata();

    for (const auto& val : _history)
      XML << "\n" << val.first << " " << val.second;
    
    XML << magnet::xml::endtag("PolarNematic");
  }
}

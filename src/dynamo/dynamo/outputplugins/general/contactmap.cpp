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

#include <dynamo/outputplugins/general/contactmap.hpp>
#include <dynamo/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPContactMap::OPContactMap(const dynamo::Simulation* t1, 
				   const magnet::xml::Node& XML):
    OutputPlugin(t1,"CollDistCheck"),
    binwidth(0.01)
  { operator<<(XML); }

  void 
  OPContactMap::initialise() 
  {}

  void 
  OPContactMap::operator<<(const magnet::xml::Node& XML)
  {}

  void 
  OPContactMap::eventUpdate(const IntEvent& eevent, 
			    const PairEventData& PDat)
  {
  }

  void 
  OPContactMap::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("ContactMap");
    XML << magnet::xml::endtag("ContactMap");
  }
}

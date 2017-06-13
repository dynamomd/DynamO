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

#include <dynamo/outputplugins/tickerproperty/craig.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <magnet/xmlwriter.hpp>
#include <fstream>
#include <sstream>

namespace dynamo {
  OPCraig::OPCraig(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"Craig")
  {
    operator<<(XML);
  }

  void 
  OPCraig::operator<<(const magnet::xml::Node& XML)
  {
    //Here is where you process options for the output plugin
  }

  void 
  OPCraig::initialise()
  {
    //This is called once, after the simulation is set up, just before
    //the first event is run.
    ticker();
  }

  void 
  OPCraig::ticker()
  {
    //This is called periodically, as set by the -t option of dynarun
  }

  void 
  OPCraig::output(magnet::xml::XmlStream& xml)
  {
    //This is your chance to output into the output.xml.bz2 file (xml arg).
    //This may be called many times if the snapshot mode is on. 
  }
}

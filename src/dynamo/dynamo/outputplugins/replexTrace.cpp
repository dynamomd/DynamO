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

#include <dynamo/outputplugins/replexTrace.hpp>
#include <dynamo/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <sstream>
#include <array>

namespace dynamo {
  OPReplexTrace::OPReplexTrace(const dynamo::Simulation* t1, const magnet::xml::Node&):
    OutputPlugin(t1,"ReplexTrace")
  {}

  void 
  OPReplexTrace::replicaExchange(OutputPlugin& OPP)
  {
    auto other = static_cast<OPReplexTrace&>(OPP);
    addPoint();
    other.addPoint();
    std::swap(Sim, other.Sim);
    addPoint();
    other.addPoint();
    std::swap(entries, other.entries);
  }

  void 
  OPReplexTrace::addPoint()
  {
    const std::array<double,3>& 
      ensembleVals(Sim->ensemble->getReducedEnsembleVals());

    std::ostringstream op;

    op << Sim->replexExchangeNumber << " "
       << Sim->systemTime / Sim->units.unitTime() << " "
       << ensembleVals[0] << " " 
       << ensembleVals[1] << " " 
       << ensembleVals[2] << "\n";

    entries.push_back(op.str());
  }

  void 
  OPReplexTrace::output(magnet::xml::XmlStream& XML)
  {
    addPoint();

    XML << magnet::xml::tag("ReplexTrace")
	<< magnet::xml::chardata();
  
    for (const std::string& str : entries)
      XML << str;

    XML << magnet::xml::endtag("ReplexTrace");

    entries.pop_back();
  }
}

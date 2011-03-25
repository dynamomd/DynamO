/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "replexTrace.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../extcode/xmlwriter.hpp"
#include <sstream>

OPReplexTrace::OPReplexTrace(const DYNAMO::SimData* t1, const XMLNode&):
  OutputPlugin(t1,"ReplexTrace")
{}

void 
OPReplexTrace::changeSystem(OutputPlugin* OPP)
{
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<OPReplexTrace*>(OPP) == NULL)
    M_throw() << "Not the correct plugin to change System with";
#endif

  addPoint();
  static_cast<OPReplexTrace*>(OPP)->addPoint();

  std::swap(Sim, static_cast<OPReplexTrace*>(OPP)->Sim);
  
  addPoint(); 
  static_cast<OPReplexTrace*>(OPP)->addPoint();

  std::swap(entries, static_cast<OPReplexTrace*>(OPP)->entries);
}

void 
OPReplexTrace::addPoint()
{
  const boost::array<double,3>& ensembleVals(Sim->ensemble->getReducedEnsembleVals());

  std::ostringstream op;

  op << Sim->dSysTime / Sim->dynamics.units().unitTime() 
	  << " " << Sim->replexExchangeNumber << " "
	  << ensembleVals[0] << " " << ensembleVals[1] << " " << ensembleVals[2] << "\n";

  entries.push_back(op.str());
}

void 
OPReplexTrace::output(xml::XmlStream& XML)
{
  addPoint();

  XML << xml::tag("ReplexTrace")
      << xml::chardata();
  
  BOOST_FOREACH(const std::string& str, entries)
    XML << str;

  XML << xml::endtag("ReplexTrace");
}

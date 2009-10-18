/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "boundedQstats.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../../schedulers/sorters/boundedPQ.hpp"

OPBoundedQStats::OPBoundedQStats(const DYNAMO::SimData* tmp, const XMLNode&):
  OPTicker(tmp,"BoundedPQstats"),
  treeSize(1)
{}

void 
OPBoundedQStats::initialise()
{  
  if (dynamic_cast<const CSSBoundedPQ*>
      (Sim->ptrScheduler->getSorter().get_ptr()) == NULL)
    D_throw() << "Not a bounded queue sorter!";

}

void
OPBoundedQStats::ticker()
{
  const CSSBoundedPQ& sorter(dynamic_cast<const CSSBoundedPQ&>
		       (*(Sim->ptrScheduler->getSorter())));
 
  treeSize.addVal(sorter.treeSize());

}

void 
OPBoundedQStats::output(xmlw::XmlStream& XML)
{
  const CSSBoundedPQ& sorter(dynamic_cast<const CSSBoundedPQ&>
			     (*(Sim->ptrScheduler->getSorter())));

  XML << xmlw::tag("boundedQstats") 
      << xmlw::attr("ExceptionEvents") << sorter.exceptionEvents()
      << xmlw::tag("CBTSize");

  treeSize.outputHistogram(XML,1.0);

  XML << xmlw::endtag("CBTSize")
    << xmlw::tag("treedist")
      << xmlw::chardata();

  if (!Sim->lNColl)
    {
      I_cerr() << "Cannot print the tree as the queue is\n"
	       << "not initialised until an event is run (i.e. N_event != 0).\n"
	       << "Continuing without tree output.";
    }
  else
    {
      const std::vector<size_t>& eventdist(sorter.getEventCounts());
      
      for (size_t i = 0; i < eventdist.size(); ++i)
	XML << i << " " 
	    << eventdist[i]
	    << "\n";
    }
  
  XML << xmlw::endtag("treedist")
      << xmlw::endtag("boundedQstats");
}

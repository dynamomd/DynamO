/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

COPBoundedQStats::COPBoundedQStats(const DYNAMO::SimData* tmp, const XMLNode&):
  COPCollTicker(tmp,"BoundedPQstats"),
  treeSize(1),
  counter(0)
{}

void 
COPBoundedQStats::initialise()
{  
  if (dynamic_cast<const CSSBoundedPQ*>
      (Sim->ptrScheduler->getSorter().get_ptr()) == NULL)
    D_throw() << "Not a bounded queue sorter!";

  eventdist.resize(dynamic_cast<const CSSBoundedPQ&>
      (*Sim->ptrScheduler->getSorter()).NLists() - 1,0);
}

void 
COPBoundedQStats::ticker()
{
  const CSSBoundedPQ& sorter(dynamic_cast<const CSSBoundedPQ&>
		       (*(Sim->ptrScheduler->getSorter())));
 
  treeSize.addVal(sorter.treeSize());

  if (!(Sim->lNColl % 100))
    {
      ++counter;
      std::vector<size_t> tmpList = sorter.getEventCounts();
      for (size_t i = 0; i < tmpList.size(); ++i)
	eventdist[i] += tmpList[i];
    }
}

void 
COPBoundedQStats::output(xmlw::XmlStream& XML)
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

  for (size_t i = 0; i < eventdist.size(); ++i)
    XML << i << " " 
	<< ((Iflt) eventdist[i])/ ((Iflt) counter)
	<< "\n";

  XML << xmlw::endtag("treedist")
      << xmlw::endtag("boundedQstats");
}

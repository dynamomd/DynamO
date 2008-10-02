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

#include "error.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../dynamics/interactions/intEvent.hpp"

COPError::COPError(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"Error"),
  allReverseColls(0), 
  strongReverseColls(0)
{}

COPError::~COPError()
{
  //std::cout << "COPError: unloaded\n";
}

void 
COPError::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  //Error Checking
  if (collision.getdt() < -0.0)
    {
      allReverseColls++;  
      if (collision.getdt() < -1e-6)
	strongReverseColls++;
    }
}

void
COPError::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Errors")
      << xmlw::tag("ReverseTimeCollisions")
      << xmlw::attr("val") << allReverseColls
      << xmlw::endtag("ReverseTimeCollisions")
      << xmlw::tag("severeTimeCollisions")
      << xmlw::attr("val") << strongReverseColls
      << xmlw::endtag("severeTimeCollisions")
      << xmlw::endtag("Errors");
}

void
COPError::periodicOutput()
{}

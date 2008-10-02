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

#include "deltaKE.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../dynamics/include.hpp"
#include "../base/is_simdata.hpp"

COPdeltaKE::COPdeltaKE(DYNAMO::SimData* tmp):
  COutputPlugin(tmp, "DeltaKE"),
  deltaKE(0.001*Sim->Dynamics.units().unitEnergy())
{}

COPdeltaKE::~COPdeltaKE()
{}

void 
COPdeltaKE::collisionUpdate(const CIntEvent &collision, 
			    const CIntEventData &preColl)
{
  deltaKE.addVal(preColl.getDeltae());
}

void
COPdeltaKE::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("DeltaKE");
  deltaKE.outputHistogram(XML, 1.0/Sim->Dynamics.units().unitEnergy());
  XML << xmlw::endtag("DeltaKE");
}

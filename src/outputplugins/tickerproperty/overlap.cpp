/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "overlap.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

OPOverlapTest::OPOverlapTest(const DYNAMO::SimData* tmp, const XMLNode&):
  OPTicker(tmp,"OverlapTester")
{}


void 
OPOverlapTest::initialise()
{
  I_cout() << "Testing for overlaps in starting configuration";
  ticker();
}

void
OPOverlapTest::output(xmlw::XmlStream&)
{
  I_cout() << "Testing for overlaps in output configuration";
  ticker();
}

void 
OPOverlapTest::ticker()
{
  for (std::vector<CParticle>::const_iterator iPtr = Sim->vParticleList.begin();
       iPtr != Sim->vParticleList.end(); ++iPtr)
    for (std::vector<CParticle>::const_iterator jPtr = iPtr + 1;
	 jPtr != Sim->vParticleList.end(); ++jPtr)
      Sim->dynamics.getInteraction(*iPtr, *jPtr)->checkOverlaps(*iPtr, *jPtr);
}

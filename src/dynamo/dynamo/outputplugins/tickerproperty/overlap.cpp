/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/outputplugins/tickerproperty/overlap.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPOverlapTest::OPOverlapTest(const dynamo::SimData* tmp, 
			       const magnet::xml::Node&):
    OPTicker(tmp,"OverlapTester")
  {}


  void 
  OPOverlapTest::initialise()
  {
    dout << "Testing for overlaps in starting configuration" << std::endl;
    ticker();
  }

  void
  OPOverlapTest::output(magnet::xml::XmlStream&)
  {
    dout << "Testing for overlaps in output configuration" << std::endl;
    ticker();
  }

  void 
  OPOverlapTest::ticker()
  {
    for (std::vector<Particle>::const_iterator iPtr = Sim->particleList.begin();
	 iPtr != Sim->particleList.end(); ++iPtr)
      for (std::vector<Particle>::const_iterator jPtr = iPtr + 1;
	   jPtr != Sim->particleList.end(); ++jPtr)
	Sim->dynamics.getInteraction(*iPtr, *jPtr)->checkOverlaps(*iPtr, *jPtr);
  }
}

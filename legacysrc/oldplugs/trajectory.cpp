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

#include "trajectory.hpp"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
namespace io = boost::iostreams;
#include <cmath>
#include "../dynamics/interactions/intEvent.hpp"
#include "../simulation/particle.hpp"
#include "../base/is_exception.hpp"

struct mempiece
{
  mempiece(Iflt ndt, unsigned long np1, unsigned long np2,
	   EEventType nType):
    dt(ndt),p1(np1),p2(np2),Type(nType) {};
  
  mempiece():
    dt(HUGE_VAL),p1(0),p2(0),Type(NONE) {};
  
  Iflt dt;
  unsigned long p1, p2;
  EEventType Type;
};

std::list<mempiece> COPTrajectory::collHistory;

void 
COPTrajectory::setTotalCollCount(unsigned long a)
{ 
  I_cout() << "Grabbing memory for trajectory";
  I_cout() << "Need " << sizeof(mempiece)*a/1048576 << "MB";
  std::cout.flush();
  collHistory.resize(a); 
  I_cout() << "Memory obtained";
}


COPTrajectory::COPTrajectory(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"Trajectory"),
  fileName("trajectory.trj"),
  simType(0)
{
  currentPos = collHistory.begin();
}

COPTrajectory::~COPTrajectory()
{}

void
COPTrajectory::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  currentPos++;
  *currentPos = mempiece(collision.getdt(), collision.getParticle1().getID(),
			 collision.getParticle2().getID(), collision.getType());
  /*  if (collision.getParticle1().getID() > collision.getParticle2().getID())
    std::cerr << collision.getParticle2().getID() << " "
	      << collision.getParticle1().getID() << "\n";
  else
    std::cerr << collision.getParticle1().getID() << " "
    << collision.getParticle2().getID() << "\n";*/
}

void
COPTrajectory::output(xmlw::XmlStream &XML)
{
  if (simType == 0)
    D_throw() << "The simulation type is not set";

  io::filtering_ostream file;
  file.push(io::bzip2_compressor());
  file.push(io::file_sink(fileName,std::ios::out|std::ios::binary|std::ios::trunc));
  
  file.write(((char *) &simType), sizeof(unsigned int));

  for (std::list<mempiece>::const_iterator iPtr = collHistory.begin();
       iPtr != currentPos; iPtr++)
    file.write(((char *) &(*iPtr)),sizeof(mempiece));
}

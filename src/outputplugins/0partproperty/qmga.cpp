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

#include "qmga.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

OPQMGA::OPQMGA(const DYNAMO::SimData* tmp, const XMLNode&):
  OPCollTicker(tmp,"OPQMGA"),
  frameCount(0)
{}

void 
OPQMGA::ticker()
{
  if (!(Sim->eventCount % 1000))
    printImage();
}

void
OPQMGA::printImage()
{
  char *fileName;

  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;

  Sim->dynamics.getLiouvillean().updateAllParticles();

  if ( asprintf(&fileName, "cnf.%04d", frameCount++) < 0)
    D_throw() << "asprintf error in QMGA";
  
  std::ofstream of(fileName);
  
  if (!of.is_open())
    D_throw() << "Could not open QMGA file for writing";

  of << Sim->N << "\n"
     << Sim->aspectRatio[0] / Sim->dynamics.units().unitLength() << "\n"
     << Sim->aspectRatio[1] / Sim->dynamics.units().unitLength() << "\n"
     << Sim->aspectRatio[2] / Sim->dynamics.units().unitLength() << "\n"
     << "0.0 0.0\n";
  
  std::list<unsigned long> tmpList;

  unsigned int i = 0;
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      Vector  pos = part.getPosition();
      Sim->dynamics.BCs().applyBC(pos);
      pos /= Sim->dynamics.units().unitLength();
      
      of << pos[0] << " " << pos[1] << " " << pos[2]
	 << "0 0 0 0.0 1.0 0.0 0 0 0 " << part.getID()
	 << " " << i << "\n";
    }
  
  of.close();
}

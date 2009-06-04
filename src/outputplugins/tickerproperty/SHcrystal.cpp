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

#include "SHcrystal.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"

COPSHCrystal::COPSHCrystal(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"SHCrystal")
{}

void 
COPSHCrystal::initialise() 
{ 
  double smallestlength = HUGE_VAL;
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& pGlob, Sim->Dynamics.getGlobals())
    ticker(); 
}
void 
COPSHCrystal::ticker()
{
  
}

std::complex<float> 
localq(const CParticle& part, int l, int m)
{
  
}

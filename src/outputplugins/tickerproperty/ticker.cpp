
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

#include "ticker.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/systems/sysTicker.hpp"

COPTicker::COPTicker(const DYNAMO::SimData* t1,const char *t2):
  COutputPlugin(t1,t2)
{}

Iflt 
COPTicker::getTickerTime() const
{
  try {
    return dynamic_cast<const CSTicker&>(*Sim->Dynamics.getSystem("SystemTicker")).getPeriod();
  } catch (const std::bad_cast&)
    {
      D_throw() << "Could not upcast the SystemTicker system event to CSTicker, have you named a system as SystemTicker?";
    }
}

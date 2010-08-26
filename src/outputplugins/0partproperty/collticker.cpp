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

#include "collticker.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"

OPCollTicker::OPCollTicker(const DYNAMO::SimData* t1,const char *t2, unsigned char order):
  OutputPlugin(t1,t2,order)
{}

void 
OPCollTicker::eventUpdate(const IntEvent &event, 
		    const PairEventData &) 
{
  stream(event.getdt());
  ticker();
}

void 
OPCollTicker::eventUpdate(const GlobalEvent &event, const NEventData&) 
{
  stream(event.getdt());
  ticker();
}

void 
OPCollTicker::eventUpdate(const LocalEvent &event, const NEventData&) 
{
  stream(event.getdt());
  ticker();
}

void 
OPCollTicker::eventUpdate(const CSystem&, const NEventData&, const Iflt& dt)
{
  stream(dt);
  ticker();
}

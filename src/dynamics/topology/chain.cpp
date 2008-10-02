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

#include "chain.hpp"
#include "../../extcode/xmlwriter.hpp"
#include <boost/foreach.hpp>

CTChain::CTChain(const XMLNode& XML, DYNAMO::SimData* Sim, unsigned int ID):
  CTopology(Sim, ID)
{
  CTopology::operator<<(XML);
  
  size_t Clength = (*ranges.begin())->size();
  BOOST_FOREACH(const smrtPlugPtr<CRange>& nRange, ranges)
    if (nRange->size() != Clength)
      I_throw() << "Size mismatch in loading one of the ranges in Chain topology \"" 
		<< spName << "\"";
}

CTChain::CTChain(DYNAMO::SimData* Sim, unsigned int ID, std::string nName):
  CTopology(Sim,ID)
{
  spName = nName;
}

void 
CTChain::outputXML(xmlw::XmlStream& XML) const 
{
  XML << xmlw::attr("Type") << "Chain";
  CTopology::outputXML(XML);
}


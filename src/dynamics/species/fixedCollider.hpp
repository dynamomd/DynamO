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

#pragma once

#include "../../extcode/xmlwriter.hpp"

class SpFixedCollider:public Species
{
public:  
  SpFixedCollider(DYNAMO::SimData* sim, std::string name, 
		  CRange* r, std::string nName, 
		  unsigned int ID, std::string nIName="Bulk"):
    Species(sim, name, r, 0.0, nName, ID, nIName)
  {}
  
  SpFixedCollider(const XMLNode& XML, DYNAMO::SimData* Sim, unsigned int nID):
    Species(XML, Sim, nID)
  { operator<<(XML); }
  

  virtual void operator<<(const XMLNode& XML)
  { Species::operator<<(XML); }

  virtual SpFixedCollider* Clone() const { return new SpFixedCollider(*this); }

protected:

  virtual void outputXML(xml::XmlStream& XML) const
  {
    XML << xml::attr("Name") << spName
	<< xml::attr("IntName") << intName
	<< xml::attr("Type") << "FixedCollider"
	<< range;
  }
};

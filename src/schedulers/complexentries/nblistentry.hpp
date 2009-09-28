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

#ifndef CSCENBList_H
#define CSCENBList_H

#include "entry.hpp"

class CSCENBList: public CSCEntry
{
public:
  CSCENBList(const XMLNode&, DYNAMO::SimData* const);
  
  virtual void initialise();

  virtual void operator<<(const XMLNode&);

  virtual void getParticleNeighbourhood(const CParticle&, 
					const CGNeighbourList::nbHoodFunc&) const;

  virtual void getParticleLocalNeighbourhood(const CParticle&, 
					     const CGNeighbourList::nbHoodFunc&) const;

  virtual CSCEntry* Clone() const { return new CSCENBList(*this); }

protected:

  virtual void outputXML(xmlw::XmlStream&) const;
  
  std::string name;
  size_t nblistID;
};

#endif

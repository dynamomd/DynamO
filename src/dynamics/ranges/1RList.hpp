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

#ifndef C1RList_H
#define C1RList_H

#include "1range.hpp"
#include <vector>

class CRList: public CRange
{
public:
  CRList(const XMLNode&);

  CRList() {}

  virtual CRange* Clone() const { return new CRList(*this); };

  virtual bool isInRange(const Particle &) const;

  //The data output classes
  virtual void operator<<(const XMLNode&);
  
  virtual unsigned long size() const { return IDs.size(); };

  virtual iterator begin() const { return CRange::iterator(0, this); }

  virtual iterator end() const { return CRange::iterator(IDs.size(), this); }

  virtual unsigned long operator[](unsigned long i) const { return IDs[i]; }

  virtual unsigned long at(unsigned long i) const { return IDs.at(i); }

protected:
  virtual const unsigned long& getIteratorID(const unsigned long &i) const { return IDs[i]; }

  virtual void outputXML(xml::XmlStream&) const;

  std::vector<unsigned long> IDs;
};

#endif


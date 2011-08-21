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

#pragma once
#include <dynamo/dynamics/ranges/2range.hpp>
#include <map>
#include <list>

namespace dynamo {
  class C2RList:public C2Range
  {
  public:
    C2RList(const magnet::xml::Node&);
    C2RList() {}

    virtual C2Range* Clone() const 
    { return new C2RList(*this); };

    virtual bool isInRange(const Particle&, const Particle&) const;

    void addPair(unsigned long, unsigned long);
  
    virtual void operator<<(const magnet::xml::Node&);

    const std::map<unsigned long, std::list<unsigned long> >& getPairMap() const;
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    std::map<unsigned long, std::list<unsigned long> > pairmap;
  };
}

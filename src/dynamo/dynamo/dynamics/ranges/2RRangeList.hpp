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
#include <dynamo/base.hpp>
#include <tr1/memory>
#include <list>

namespace dynamo {
  class C2RRangeList:public C2Range, dynamo::SimBase_const
  {
  public:
    C2RRangeList(const magnet::xml::Node&, const dynamo::SimData*);
    C2RRangeList(const dynamo::SimData* nSim):SimBase_const(nSim,"C2RRangeList") {}

    virtual bool isInRange(const Particle&, const Particle&) const;

    void addRange(C2Range* nRange)
    { ranges.push_back(std::tr1::shared_ptr<C2Range>(nRange)); }
  
    virtual void operator<<(const magnet::xml::Node&);

    const std::list<std::tr1::shared_ptr<C2Range> >& getRanges() const;
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    std::list<std::tr1::shared_ptr<C2Range> > ranges;
  };
}

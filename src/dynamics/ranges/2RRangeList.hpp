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

#ifndef C2RRangeList_H
#define C2RRangeList_H

#include <list>
#include "2range.hpp"
#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"

class C2RRangeList:public C2Range, DYNAMO::SimBase_const
{
public:
  C2RRangeList(const XMLNode&, const DYNAMO::SimData*);
  C2RRangeList(const DYNAMO::SimData* nSim):SimBase_const(nSim,"C2RRangeList",IC_red) {}

  virtual C2Range* Clone() const 
  { return new C2RRangeList(*this); };

  virtual bool isInRange(const CParticle&, const CParticle&) const;

  void addRange(C2Range* nRange)
  { ranges.push_back(smrtPlugPtr<C2Range>(nRange)); }
  
  virtual void operator<<(const XMLNode&);

  const std::list<smrtPlugPtr<C2Range> >& getRanges() const;
  
protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  std::list<smrtPlugPtr<C2Range> > ranges;
};

#endif

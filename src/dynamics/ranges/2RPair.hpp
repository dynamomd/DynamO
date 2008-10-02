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

#ifndef C2RPair_H
#define C2RPair_H

#include "1range.hpp"
#include "2range.hpp"
#include "../../datatypes/pluginpointer.hpp"

class C2RPair:public C2Range
{
public:
  C2RPair(const XMLNode&, const DYNAMO::SimData*);

  C2RPair(CRange* r1, CRange* r2 ):range1(r1),range2(r2) {}
  
  virtual C2Range* Clone() const 
  { return new C2RPair(*this); };

  virtual bool isInRange(const CParticle&, const CParticle&) const;
  
  virtual void operator<<(const XMLNode&);
  
protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  smrtPlugPtr<CRange> range1;
  smrtPlugPtr<CRange> range2;
};

#endif

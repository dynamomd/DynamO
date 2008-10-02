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

#ifndef C2RSingle_H
#define C2RSingle_H

#include "1range.hpp"
#include "2range.hpp"
#include "../../datatypes/pluginpointer.hpp"

class C2RSingle:public C2Range
{
public:
  C2RSingle(const XMLNode&, const DYNAMO::SimData*);

  C2RSingle(CRange* r1):range(r1) {}
  
  virtual C2Range* Clone() const 
  { return new C2RSingle(*this); };

  virtual bool isInRange(const CParticle&, const CParticle&) const;
  
  virtual void operator<<(const XMLNode&);
  
  const smrtPlugPtr<CRange>& getRange() const { return range; }

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  smrtPlugPtr<CRange> range;
};

#endif

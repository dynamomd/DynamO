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
#include "1range.hpp"
#include "2range.hpp"
#include <magnet/cloneptr.hpp>

class C2RPair:public C2Range
{
public:
  C2RPair(const magnet::xml::Node&, const dynamo::SimData*);

  C2RPair(CRange* r1, CRange* r2 ):range1(r1),range2(r2) {}
  
  virtual C2Range* Clone() const 
  { return new C2RPair(*this); };

  virtual bool isInRange(const Particle&, const Particle&) const;
  
  virtual void operator<<(const magnet::xml::Node&);
  
protected:
  virtual void outputXML(magnet::xml::XmlStream&) const;

  magnet::ClonePtr<CRange> range1;
  magnet::ClonePtr<CRange> range2;
};

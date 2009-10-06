/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef SW_Units_H
#define SW_Units_H

#include "units.hpp"

class CUSW: public CUnits
{
 public:
  CUSW(const DYNAMO::SimData*); 

  CUSW(Iflt, Iflt, const DYNAMO::SimData*);

  CUSW(const XMLNode&, const DYNAMO::SimData*);

  virtual ~CUSW();

  virtual Iflt unitLength() const;

  virtual void setUnitLength(Iflt);

  virtual Iflt unitTime() const;
  
  virtual CUnits* Clone() const;
  
  virtual void operator<<(const XMLNode&);

  virtual void rescaleLength(Iflt);

 protected:
  virtual void outputXML(xmlw::XmlStream &) const;

  Iflt UnitOfEnergy;
  Iflt UnitOfLength;
};

#endif

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

#ifndef Elas_Units_H
#define Elas_Units_H

#include "units.hpp"

class CUElastic: public CUnits
{
 public:
  CUElastic(const DYNAMO::SimData*);
  
  CUElastic(Iflt, const DYNAMO::SimData*);

  CUElastic(const XMLNode&, const DYNAMO::SimData*);

  virtual ~CUElastic();

  virtual Iflt unitLength() const;

  virtual void setUnitLength(Iflt);

  virtual Iflt unitTime() const;

  virtual void rescaleLength(Iflt);
  
  virtual CUnits* Clone() const;
  
  virtual void operator<<(const XMLNode&);

 protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  Iflt UnitOfLength;
};

#endif

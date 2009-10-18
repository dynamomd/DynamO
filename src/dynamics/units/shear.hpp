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

#ifndef Shear_Units_H
#define Shear_Units_H

#include "elastic.hpp"

class CUShear: public CUElastic
{
public:
  CUShear(const DYNAMO::SimData* tmp): 
    CUElastic(1.0, tmp)
  {
    I_cout() << "Elastic units loaded";
  }
  
  CUShear(const Iflt& length, const DYNAMO::SimData* tmp):
    CUElastic(length, tmp)
  {
    I_cout() << "Elastic units loaded";
  }
  
  CUShear(const XMLNode& XML, const DYNAMO::SimData* tmp):
    CUElastic(1.0, tmp)
  { 
    operator<<(XML); 
    I_cout() << "Elastic units loaded";
  }
  
  virtual ~CUShear() {}
  
  virtual Iflt unitTime() const
  { return 1.0 / ShearRate; }
  
  virtual Units* Clone() const
  { return new CUShear(*this); }
  
  void operator<<(const XMLNode&);

private:
  void outputXML(xmlw::XmlStream&) const;
};

#endif

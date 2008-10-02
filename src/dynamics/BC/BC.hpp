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

#ifndef BC_H
#define BC_H

#include "../../datatypes/vector.hpp"

class XMLNode;
namespace xmlw
{
class XmlStream;
}
class CIntEvent;
namespace DYNAMO
{
  class SimData;
}

class CBC
{
 public:
  virtual ~CBC() {};

  virtual void setPBC(CVector<> &)const = 0;

  virtual void setPBC(CVector<> &, CVector<> &) const = 0;

  virtual void setPBC(CVector<> &posVec, Iflt) const 
    { setPBC(posVec); }

  virtual void update(const Iflt&) {};

  virtual void operator<<(const XMLNode &) = 0;

  virtual CBC* Clone () const = 0;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CBC&);

  static CBC* loadClass(const XMLNode& ,DYNAMO::SimData*);

 protected:
  virtual void outputXML(xmlw::XmlStream &XML) const = 0;

  virtual void rounding(CVector<> &) const = 0;
  
  DYNAMO::SimData* Sim;
};

#endif

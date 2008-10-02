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

#ifndef PBC_H
#define PBC_H

#include "shapes.hpp"
#include "../../base/is_base.hpp"

class CSPBC: public CSqBC, public DYNAMO::Base_Class
{
public:
  CSPBC(DYNAMO::SimData*);

  void setPBC(CVector<>& pos) const
  { rounding(pos); }

  void setPBC(CVector<>& pos, CVector<>&) const
  { rounding(pos); }

  inline virtual void outputXML(xmlw::XmlStream &) const;  
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

class CRPBC: public CRectBC, public DYNAMO::Base_Class
{
public:
  CRPBC(DYNAMO::SimData*);

  void setPBC(CVector<>& pos) const
  { rounding(pos); }
  
  void setPBC(CVector<>& pos, CVector<>&) const
  { rounding(pos); }

  virtual void outputXML(xmlw::XmlStream&) const;
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

class CRNoXPBC: public CRectBC, public DYNAMO::Base_Class
{
public:
  CRNoXPBC(DYNAMO::SimData*);

  void setPBC(CVector<>& pos) const
  { 
    Iflt x = pos[0];
    rounding(pos); 
    pos[0] = x;
  }
  
  void setPBC(CVector<>& pos, CVector<>&) const
  { 
    Iflt x = pos[0];
    rounding(pos); 
    pos[0] = x;
  }

  virtual void outputXML(xmlw::XmlStream&) const;
  virtual void operator<<(const XMLNode&);
  virtual CBC* Clone () const;
};

#endif

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

#ifndef LEBC_H
#define LEBC_H

#include "shapes.hpp"
#include "../../base/is_base.hpp"

class CRLEBC: virtual public CRectBC, public DYNAMO::Base_Class
{
 public:
  CRLEBC(DYNAMO::SimData*);

  CRLEBC(const XMLNode&, DYNAMO::SimData*);

  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void operator<<(const XMLNode&);

  virtual CBC* Clone () const;

  virtual void setPBC(CVector<>&) const; 

  virtual void setPBC(CVector<>&, CVector<>&) const;

  virtual void setPBC(CVector<>&, Iflt dt) const;

  virtual void update(const Iflt&);

 protected:  
  Iflt dxd;
};

class CSLEBC: virtual public CSqBC, public DYNAMO::Base_Class
{
 public:
  CSLEBC(DYNAMO::SimData*);

  CSLEBC(const XMLNode&, DYNAMO::SimData*);

  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void operator<<(const XMLNode&);

  virtual CBC* Clone () const;

  virtual void setPBC(CVector<>&) const; 

  virtual void setPBC(CVector<>&, CVector<>&) const;

  virtual void setPBC(CVector<>&, Iflt dt) const;

  virtual void update(const Iflt&);

 protected:  
  Iflt dxd;
};

#endif

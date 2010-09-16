/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSysNull_HPP
#define CSysNull_HPP

#include "system.hpp"
#include "../../base/is_exception.hpp"
#include "../NparticleEventData.hpp"

class CSysNull: public System
{
public:
  CSysNull(DYNAMO::SimData* tmp): System(tmp) 
  {  sysName = "NULL"; }
  
  virtual System* Clone() const { return new CSysNull(*this); }

  virtual void stream(Iflt) {}

  virtual NEventData runEvent()
  { M_throw() << "You're running the null system event"; }

  virtual void initialise() {}

  virtual void operator<<(XMLNode&) {}

  virtual bool operator<(const IntEvent&) const { return false; }

  virtual bool operator<(const GlobalEvent&) const { return false; }

  virtual bool operator<(const System&) const { return false; }

protected:
  virtual void outputXML(xml::XmlStream&) const {}
};

#endif

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

#ifndef CLWall_HPP
#define CLWall_HPP

#include "local.hpp"

class CLWall: public CLocal
{
public:
  CLWall(const XMLNode&, DYNAMO::SimData*);
  CLWall(DYNAMO::SimData*, Iflt, CVector<>, CVector<>, 
	 std::string, CRange*);

  virtual ~CLWall() {}

  virtual CLocal* Clone() const { return new CLWall(*this); };

  virtual CLocalEvent getEvent(const CParticle&) const;

  virtual void runEvent(const CParticle&, const CLocalEvent&) const;
  
  virtual bool isInCell(const CVector<>&, const CVector<>&) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  CVector<> vNorm;
  CVector<> vPosition;
  Iflt e;
};

#endif

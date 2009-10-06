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

#ifndef CLAndersenWall_HPP
#define CLAndersenWall_HPP

#include "local.hpp"
#include "../../datatypes/vector.hpp"

class CLAndersenWall: public CLocal
{
public:
  CLAndersenWall(const XMLNode&, DYNAMO::SimData*);

  CLAndersenWall(DYNAMO::SimData*, Iflt, Vector , Vector , 
		 std::string, CRange*);

  virtual ~CLAndersenWall() {}

  virtual CLocal* Clone() const { return new CLAndersenWall(*this); }

  virtual CLocalEvent getEvent(const CParticle &) const;

  virtual void runEvent(const CParticle&, const CLocalEvent&) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

  virtual bool isInCell(const Vector &, const Vector &) const;

  virtual void write_povray_info(std::ostream&) const;

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  Vector  vNorm;
  Vector  vPosition;
  Iflt sqrtT;
};

#endif

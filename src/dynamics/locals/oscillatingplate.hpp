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

#ifndef CLOscillatingPlate_HPP
#define CLOscillatingPlate_HPP

#include "local.hpp"

class CLOscillatingPlate: public CLocal
{
public:
  CLOscillatingPlate(const XMLNode&, DYNAMO::SimData*);
  CLOscillatingPlate(DYNAMO::SimData*, Vector, Vector, Iflt, 
		     Iflt, Iflt, Iflt, std::string, CRange*);

  virtual ~CLOscillatingPlate() {}

  virtual CLocal* Clone() const { return new CLOscillatingPlate(*this); };

  virtual CLocalEvent getEvent(const CParticle&) const;

  virtual void runEvent(const CParticle&, const CLocalEvent&) const;
  
  virtual bool isInCell(const Vector &, const Vector &) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

  virtual void write_povray_info(std::ostream&) const;

  Vector getPosition() const;

protected:
  virtual void outputXML(xmlw::XmlStream&) const;
  
  Vector rw0;
  Vector nhat;
  Iflt omega0;
  Iflt sigma;
  Iflt e;
  Iflt delta;
};

#endif

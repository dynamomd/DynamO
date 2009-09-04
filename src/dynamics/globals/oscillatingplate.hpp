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

#ifndef CGOscillatingPlate_HPP
#define CGOscillatingPlate_HPP

#include "global.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../../datatypes/vector.hpp"
#include "../../simulation/particle.hpp"
#include <vector>

class CGOscillatingPlate: public CGlobal
{
public:
  CGOscillatingPlate(const XMLNode&, DYNAMO::SimData*);

  CGOscillatingPlate(DYNAMO::SimData*, Iflt, Iflt, Iflt, Iflt, const std::string&);

  virtual ~CGOscillatingPlate() {}

  virtual CGlobal* Clone() const 
  { 
    return new CGOscillatingPlate(*this); 
  }

  virtual CGlobEvent getEvent(const CParticle &) const;

  virtual void runEvent(const CParticle&) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

  virtual void outputXML(xmlw::XmlStream& XML) const;

protected:

  Iflt x0;
  Iflt xi;
  Iflt omega0;
  Iflt sigma;
};

#endif

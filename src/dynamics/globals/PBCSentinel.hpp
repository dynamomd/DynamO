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

#ifndef CGPBCSentinel_HPP
#define CGPBCSentinel_HPP

#include "global.hpp"
#include <vector>

class CGPBCSentinel: public CGlobal
{
public:
  CGPBCSentinel(const XMLNode&, DYNAMO::SimData*);

  CGPBCSentinel(DYNAMO::SimData*, const std::string&);
  
  virtual ~CGPBCSentinel() {}

  virtual CGlobal* Clone() const { return new CGPBCSentinel(*this); };

  virtual CGlobEvent getEvent(const CParticle &) const;

  virtual void runEvent(const CParticle&) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  void particlesUpdated(const CNParticleData&);

  virtual void outputXML(xmlw::XmlStream&) const;

  Iflt maxintdist;

  mutable std::vector<Iflt> cachedTimes;
};

#endif

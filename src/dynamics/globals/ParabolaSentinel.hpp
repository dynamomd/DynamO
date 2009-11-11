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

#ifndef CGParabolaSentinel_HPP
#define CGParabolaSentinel_HPP

#include "global.hpp"
#include <vector>

class CGParabolaSentinel: public CGlobal
{
public:
  CGParabolaSentinel(const XMLNode&, DYNAMO::SimData*);

  CGParabolaSentinel(DYNAMO::SimData*, const std::string&);
  
  virtual ~CGParabolaSentinel() {}

  virtual CGlobal* Clone() const { return new CGParabolaSentinel(*this); };

  virtual CGlobEvent getEvent(const CParticle &) const;

  virtual void runEvent(const CParticle&) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  void particlesUpdated(const CNParticleData&);

  virtual void outputXML(xmlw::XmlStream&) const;
  
  mutable std::vector<unsigned char> passedParabola;
};

#endif

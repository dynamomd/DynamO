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

#ifndef COP1PP_HPP
#define COP1PP_HPP

#include "../outputplugin.hpp"

class COP1PP: public COutputPlugin
{
public:
  COP1PP(const DYNAMO::SimData*, const char*, unsigned char order = 100);

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

private:
  virtual void A2ParticleChange(const C2ParticleData&);

  /* This class of plugins implement these functions */
  virtual void A1ParticleChange(const C1ParticleData&) = 0;
  virtual void stream(const Iflt&) = 0;  
};

#endif

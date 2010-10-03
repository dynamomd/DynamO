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

#ifndef OP1PP_HPP
#define OP1PP_HPP

#include "../outputplugin.hpp"

class OP1PP: public OutputPlugin
{
public:
  OP1PP(const DYNAMO::SimData*, const char*, unsigned char order = 100);

  virtual void eventUpdate(const IntEvent&, const PairEventData&);

  virtual void eventUpdate(const GlobalEvent&, const NEventData&);

  virtual void eventUpdate(const LocalEvent&, const NEventData&);

  virtual void eventUpdate(const System&, const NEventData&, const double&);

private:
  virtual void A2ParticleChange(const PairEventData&);

  /* This class of plugins implement these functions */
  virtual void A1ParticleChange(const ParticleEventData&) = 0;
  virtual void stream(const double&) = 0;  
};

#endif

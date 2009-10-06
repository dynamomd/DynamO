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

#ifndef COPSCParameter_H
#define COPSCParameter_H

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "ticker.hpp"

class COPSCParameter: public COPTicker
{
 public:
  COPSCParameter(const DYNAMO::SimData*, const XMLNode&);

  virtual COutputPlugin *Clone() const
  { return new COPSCParameter(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();
  
  virtual void output(xmlw::XmlStream&);

  virtual void operator<<(const XMLNode&);

 protected:
  
  size_t maxWaveNumber;
  size_t count;
  std::vector<Iflt> runningsum;
};

#endif


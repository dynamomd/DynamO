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

#ifndef COPMFT_H
#define COPMFT_H

#include "1partproperty.hpp"
#include <vector>
#include <boost/circular_buffer.hpp>
#include "../../datatypes/histogram.hpp"

class COPMFT: public COP1PP
{
 public:
  COPMFT(const DYNAMO::SimData*, const XMLNode&);

  void A1ParticleChange(const C1ParticleData&);

  void stream(const Iflt&) {}

  void output(xmlw::XmlStream &); 

  void periodicOutput() {}

  virtual void initialise();

  virtual COutputPlugin *Clone() const { return new COPMFT(*this); }
  
  virtual void operator<<(const XMLNode&);
  
 protected:
  size_t collisionHistoryLength;
  
  Iflt binwidth;

  //! Each particles last collision times
  std::vector<boost::circular_buffer<Iflt> > lastTime;

  //! A histogram for each species
  std::vector<std::vector<C1DHistogram> > data;
};

#endif

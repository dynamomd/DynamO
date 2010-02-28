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

#ifndef OPCollisionCorrelator_HPP
#define OPCollisionCorrelator_HPP

#include "2partproperty.hpp"
#include "../../datatypes/histogram.hpp"
#include <boost/circular_buffer.hpp>
#include <vector>

class OPCollisionCorrelator: public OP2PP
{
public:
  OPCollisionCorrelator(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();

  virtual OutputPlugin* Clone() const 
  { return new OPCollisionCorrelator(*this); }

  void output(xmlw::XmlStream &XML);

  void operator<<(const XMLNode&);

private:

  virtual void A2ParticleChange(const C2ParticleData&);

  virtual void stream(const Iflt&) {}  

  C1DHistogram freetimehist;
  std::vector<std::vector<Iflt> > lastColl;
};

#endif

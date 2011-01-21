/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef OPChatteringCorrelator_HPP
#define OPChatteringCorrelator_HPP

#include "2partproperty.hpp"
#include "../../datatypes/histogram.hpp"
#include <boost/circular_buffer.hpp>
#include <vector>

class OPChatteringCorrelator: public OP2PP
{
public:
  OPChatteringCorrelator(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();

  virtual OutputPlugin* Clone() const
  { return new OPChatteringCorrelator(*this); }

  void output(xml::XmlStream &XML);

private:

  virtual void A2ParticleChange(const PairEventData&);

  virtual void stream(const double&) {}

  C1DWeightHistogram hist;
  std::vector<std::pair<double,double> > chatterTracker;
};

#endif

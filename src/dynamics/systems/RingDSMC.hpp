/*  dynamo:- Event driven molecular dynamics simulator 
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

#pragma once
#include "system.hpp"
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <magnet/cloneptr.hpp>

class CSRingDSMC: public System
{
public:
  CSRingDSMC(const magnet::xml::Node& XML, dynamo::SimData*);

  CSRingDSMC(dynamo::SimData*, double, double, double, double, double, std::string, CRange*);
  
  virtual System* Clone() const { return new CSRingDSMC(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node&);

protected:
  virtual void outputXML(xml::XmlStream&) const;

  mutable boost::variate_generator<dynamo::baseRNG&, 
				   boost::uniform_real<> > uniformRand;  

  double tstep;
  double chi12, chi13;
  double d2;
  double diameter;
  mutable double maxprob12;
  mutable double maxprob13;
  double e;
  double factor12;
  double factor13;

  mutable unsigned long n12;
  mutable unsigned long n13;

  magnet::ClonePtr<CRange> range1;
};

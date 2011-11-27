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
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>

namespace dynamo {
  class SysRingDSMC: public System
  {
  public:
    SysRingDSMC(const magnet::xml::Node& XML, dynamo::SimData*);

    SysRingDSMC(dynamo::SimData*, double, double, double, double, double, std::string, Range*);
  
    virtual void runEvent() const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

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

    shared_ptr<Range> range1;
  };
}

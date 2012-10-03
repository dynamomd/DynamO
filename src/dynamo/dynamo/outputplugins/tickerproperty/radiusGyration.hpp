/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <magnet/math/histogram.hpp>
#include <magnet/math/vector.hpp>
#include <list>

class Range;

namespace dynamo {
  class TChain;

  class OPRGyration: public OPTicker
  {
  public:
    OPRGyration(const dynamo::Simulation*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();

    virtual void changeSystem(OutputPlugin*);

    virtual void output(magnet::xml::XmlStream&);

    struct molGyrationDat
    {
      Vector  EigenVal;
      Vector  EigenVec[3];
      Vector  MassCentre;
    };
  
    static molGyrationDat getGyrationEigenSystem(const shared_ptr<Range>&, const dynamo::Simulation*);

    static Vector  NematicOrderParameter(const std::list<Vector  >&);

    virtual void operator<<(const magnet::xml::Node&);
  
  protected:

    struct CTCdata
    {
      const TChain* chainPtr;
      std::vector<magnet::math::Histogram<> > gyrationRadii;
      std::vector<magnet::math::Histogram<> > nematicOrder;

      CTCdata(const TChain* ptr, double binwidth1, double binwidth2, double binwidth3):
	chainPtr(ptr)
      {
	for (size_t i = 0; i < NDIM; i++)
	  {
	    gyrationRadii.push_back(magnet::math::Histogram<>(binwidth1));
	    nematicOrder.push_back(magnet::math::Histogram<>(binwidth2));
	  }
      }

    };

    std::list<CTCdata> chains;

    double binwidth1, binwidth2, binwidth3;  
  };
}

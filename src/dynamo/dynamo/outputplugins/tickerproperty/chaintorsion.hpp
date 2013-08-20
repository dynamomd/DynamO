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
#include <list>

namespace dynamo {
  class TChain;

  class OPCTorsion: public OPTicker
  {
  public:
    OPCTorsion(const dynamo::Simulation*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();

    virtual void replicaExchange(OutputPlugin&);

    virtual void output(magnet::xml::XmlStream&);
  
  protected:

    struct CTCdata
    {
      const TChain* chainPtr;
      magnet::math::Histogram<> gammaMol;
      magnet::math::Histogram<> gammaSys;
      magnet::math::Histogram<> f;
      CTCdata(const TChain* ptr, double binwidth1, double binwidth2, double binwidth3):
	chainPtr(ptr), gammaMol(binwidth1),
	gammaSys(binwidth2), f(binwidth3)
      {}

    };

    std::list<CTCdata> chains;
  };
}

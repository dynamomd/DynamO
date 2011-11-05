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
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <magnet/math/histogram.hpp>
#include <boost/shared_array.hpp>

namespace dynamo {
  class CTChain;
  class CRange;

  class OPCContactMap: public OPTicker
  {
  public:
    OPCContactMap(const dynamo::SimData*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();

    virtual void changeSystem(OutputPlugin*);

    virtual void output(magnet::xml::XmlStream&);
  
  protected:

    struct Cdata
    {
      Cdata(const CTChain*, unsigned long);

      const CTChain* chainPtr;
      boost::shared_array<unsigned long> array;
      unsigned long counter;
      unsigned long chainlength;
    };

    std::list<Cdata> chains;

  };
}

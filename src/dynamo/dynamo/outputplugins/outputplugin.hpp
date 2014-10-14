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
#include <dynamo/base.hpp>
#include <dynamo/eventtypes.hpp>

namespace magnet { namespace xml { class Node; class XmlStream; } }

namespace dynamo {
  struct XMLNode;
  class PairEventData;
  class ParticleEventData;
  class NEventData;
  class System;
  class IDRange;

  class OutputPlugin: public dynamo::SimBase_const
  {
  public:
    OutputPlugin(const dynamo::Simulation*, const char*, unsigned char order=100);
  
    inline virtual ~OutputPlugin() {}
  
    virtual void initialise() = 0;
  
    virtual void eventUpdate(const Event&, const NEventData&) = 0;
  
    virtual void output(magnet::xml::XmlStream&);
  
    virtual void periodicOutput();
  
    static shared_ptr<OutputPlugin> getPlugin(const magnet::xml::Node&, const dynamo::Simulation*);
    static shared_ptr<OutputPlugin> getPlugin(const std::string, const dynamo::Simulation*);
  
    inline bool operator<(const OutputPlugin& OP) const
    { return updateOrder < OP.updateOrder; }
  
    inline bool operator>(const OutputPlugin& OP) const
    { return updateOrder > OP.updateOrder; }
  
    virtual void replicaExchange(OutputPlugin&) = 0;
  
    virtual void temperatureRescale(const double&) {}
  
  protected:
    std::ostream& I_Pcout() const;
  
    // This sets the order in which these things are updated
    // 0 is first
    // 100 is default
    // 250 is last
    //
    // Lets other plugins take data from plugins before/after they are updated
    unsigned char updateOrder;
  };
}

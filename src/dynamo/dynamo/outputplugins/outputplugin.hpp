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
#include <dynamo/base.hpp>
#include <dynamo/base/is_simdata.hpp>

namespace xml { class XmlStream; }
namespace magnet { namespace xml { class Node; } }

namespace dynamo {
  class IntEvent;
  class GlobalEvent;
  struct XMLNode;
  class PairEventData;
  class ParticleEventData;
  class NEventData;
  class System;

  class OutputPlugin: public dynamo::SimBase_const
  {
  public:
    OutputPlugin(const dynamo::SimData*, const char*, unsigned char order=100);
  
    inline virtual ~OutputPlugin() {}
  
    virtual void initialise() = 0;
  
    virtual void eventUpdate(const IntEvent&, const PairEventData&) = 0;
  
    virtual void eventUpdate(const GlobalEvent&, const NEventData&) = 0;

    virtual void eventUpdate(const LocalEvent&, const NEventData&) = 0;

    virtual void eventUpdate(const System&, const NEventData&, const double&) = 0;
  
    virtual OutputPlugin *Clone() const = 0;
  
    virtual void output(magnet::xml::XmlStream&);
  
    virtual void periodicOutput();
  
    static OutputPlugin* getPlugin(const magnet::xml::Node&, const dynamo::SimData*);
    static OutputPlugin* getPlugin(const std::string, const dynamo::SimData*);
  
    inline bool operator<(const OutputPlugin& OP) const
    { return updateOrder < OP.updateOrder; }
  
    inline bool operator>(const OutputPlugin& OP) const
    { return updateOrder > OP.updateOrder; }
  
    virtual void changeSystem(OutputPlugin*) 
    { M_throw() << "This plugin hasn't been prepared for changes of system\n Plugin " <<  name; }
  
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

  private:

    template<class T> static OutputPlugin* 
    testGeneratePlugin(const dynamo::SimData*, const magnet::xml::Node&);
  };
}

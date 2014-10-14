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

#include <dynamo/outputplugins/outputplugin.hpp>
#include <vector>
#include <string>

namespace dynamo {
  class OPReplexTrace: public OutputPlugin
  {
  public:
    OPReplexTrace(const dynamo::Simulation*, const magnet::xml::Node&);

    void eventUpdate(const Event&, const NEventData&) {}

    virtual void initialise() { addPoint(); }

    virtual void output(magnet::xml::XmlStream&);

    virtual void replicaExchange(OutputPlugin&); 

  private:
    void addPoint();

    mutable std::vector<std::string> entries;
  };
}

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

namespace dynamo {
  /*! \brief An output plugin marker class for periodically 'ticked'
   * plugins, ticked by the SysTicker class.
   *
   * This class doesn't require any Liouvillean::updateParticle or
   * Liouvillean::updateAllParticles as this is done in the SysTicker
   * class. This is optimal as most ticker plugins need it anyway
   */
  class OPTicker: public OutputPlugin
  {
  public:
    OPTicker(const dynamo::SimData*, const char*);

    //Non virtual to warn if you use them,
    void eventUpdate(const IntEvent&, const PairEventData&) {}
    void eventUpdate(const GlobalEvent&, const NEventData&) {}
    void eventUpdate(const LocalEvent&, const NEventData&) {}
    void eventUpdate(const System&, const NEventData&, const double&) {}

    virtual void output(magnet::xml::XmlStream&) {}

    virtual void ticker() = 0;
  
    virtual void periodicOutput() {}

  protected:

    double getTickerTime() const;
  };
}

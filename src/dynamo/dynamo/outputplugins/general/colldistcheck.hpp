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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/datatypes/histogram.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <map>

namespace dynamo {
  using namespace EventTypeTracking;

  class OPCollMatrix;

  class OPCollDistCheck: public OutputPlugin
  {
  public:
    OPCollDistCheck(const dynamo::SimData*, 
		    const magnet::xml::Node&);

    ~OPCollDistCheck();

    void eventUpdate(const IntEvent&, const PairEventData&);

    void eventUpdate(const GlobalEvent&, const NEventData&);

    void eventUpdate(const LocalEvent&, const NEventData&);
  
    void eventUpdate(const System&, const NEventData&, const double&);

    OutputPlugin *Clone() const { return new OPCollDistCheck(*this); }

    virtual void initialise();

    virtual void output(magnet::xml::XmlStream&);

    virtual void operator<<(const magnet::xml::Node&);

  private:

    typedef std::pair<classKey, EEventType> eventKey;
  
    std::map<eventKey, C1DHistogram> distList;
  
    double binwidth;
  };
}

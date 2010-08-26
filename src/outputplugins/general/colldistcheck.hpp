/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef OPCollDistCheck_HPP
#define OPCollDistCheck_HPP

#include "../outputplugin.hpp"
#include <map>
#include "../../datatypes/histogram.hpp"
#include "../eventtypetracking.hpp"

using namespace EventTypeTracking;

class OPCollMatrix;

class OPCollDistCheck: public OutputPlugin
{
public:
  OPCollDistCheck(const DYNAMO::SimData*, const XMLNode&);

  ~OPCollDistCheck();

  void eventUpdate(const IntEvent&, const PairEventData&);

  void eventUpdate(const GlobalEvent&, const NEventData&);

  void eventUpdate(const LocalEvent&, const NEventData&);
  
  void eventUpdate(const CSystem&, const NEventData&, const Iflt&);

  OutputPlugin *Clone() const { return new OPCollDistCheck(*this); }

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

  virtual void operator<<(const XMLNode&);

private:

  typedef std::pair<classKey, EEventType> eventKey;
  
  std::map<eventKey, C1DHistogram> distList;
  
  Iflt binwidth;
};

#endif

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
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/interactions/captures.hpp>
#include <unordered_map>

namespace dynamo {
  /*! \brief A Dynamics which implements Newtonian dynamics, but with
    a deformed energy landscape set controlled through an interaction
    contact map.
   */
  class DynNewtonianMCCMap: public DynNewtonian
  {
    struct WData{
      WData(): _distance(0), _wval(0) {}
      WData(size_t distance, double wval): _distance(distance), _wval(wval) {}
      size_t _distance;
      double _wval;
    };

    typedef std::vector<std::pair<detail::CaptureMapKey, WData> > WContainer;
    
    WContainer _W;

    std::string _interaction_name;
    std::shared_ptr<ICapture> _interaction;
  public:
    DynNewtonianMCCMap(dynamo::Simulation* tmp, const magnet::xml::Node&);

    virtual PairEventData SphereWellEvent(const IntEvent&, const double&, const double&, size_t) const;
    virtual NEventData multibdyWellEvent(const IDRange&, const IDRange&, const double&, const double&, EEventType&) const;
    virtual void initialise();
    virtual void replicaExchange(Dynamics& oDynamics);

    double W(const detail::CaptureMap& map) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream& ) const;
  };
}

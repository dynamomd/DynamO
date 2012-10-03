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
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/math/timeaveragedproperty.hpp>
#include <magnet/math/correlators.hpp>
#include <map>
#include <ctime>

namespace dynamo {
  using namespace EventTypeTracking;

  class OPMisc: public OutputPlugin
  {

  public:
    OPMisc(const dynamo::Simulation*, const magnet::xml::Node&);
  
    virtual void initialise();
  
    virtual void eventUpdate(const IntEvent&, const PairEventData&);
  
    virtual void eventUpdate(const GlobalEvent&, const NEventData&);

    virtual void eventUpdate(const LocalEvent&, const NEventData&);
  
    virtual void eventUpdate(const System&, const NEventData&, 
			     const double&);
  
    void output(magnet::xml::XmlStream &); 
  
    void periodicOutput();
  
    double getMFT() const;
  
    void changeSystem(OutputPlugin*);

    double getEventsPerSecond() const;
    double getSimTimePerSecond() const;

    void temperatureRescale(const double&);

    double getMeankT() const;
    double getMeanSqrkT() const;
    double getCurrentkT() const;

    double getMeanUConfigurational() const;
    double getMeanSqrUConfigurational() const;
    inline double getConfigurationalU() const { return _internalE.current(); }

  protected:
    void newEvent(const size_t&, const EEventType&, const classKey&);
  
    typedef std::pair<classKey, EEventType> EventKey;
    std::map<EventKey, size_t> _counters;

    void stream(double dt);
    void eventUpdate(const NEventData&);

    std::time_t tstartTime;
    timespec acc_tstartTime;

    unsigned long _dualEvents;  
    unsigned long _singleEvents;
    unsigned long _virtualEvents;    
    size_t _reverseEvents;

    magnet::math::TimeAveragedProperty<double> _KE;
    magnet::math::TimeAveragedProperty<double> _internalE;
    magnet::math::TimeAveragedProperty<Vector> _sysMomentum;
    magnet::math::TimeAveragedProperty<Matrix> _kineticP;

    magnet::math::LogarithmicTimeCorrelator<Vector> _thermalConductivity;
    magnet::math::LogarithmicTimeCorrelator<Matrix> _viscosity;
    std::vector<magnet::math::LogarithmicTimeCorrelator<Vector> > _thermalDiffusion;
    std::vector<magnet::math::LogarithmicTimeCorrelator<Vector> > _mutualDiffusion;

    std::vector<double> _speciesMasses;
    std::vector<Vector> _speciesMomenta;
    double _systemMass;

    Matrix collisionalP;
  };
}

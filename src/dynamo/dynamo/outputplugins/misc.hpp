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
#include <chrono>
#include <map>

namespace dynamo {
  using namespace EventTypeTracking;

  class OPMisc: public OutputPlugin
  {

  public:
    OPMisc(const dynamo::Simulation*, const magnet::xml::Node&);
  
    virtual void initialise();
  
    virtual void eventUpdate(const Event&, const NEventData&);
  
    void output(magnet::xml::XmlStream &); 
  
    void periodicOutput();
  
    double getMFT() const;
  
    void replicaExchange(OutputPlugin&);

    double getDuration() const;
    double getEventsPerSecond() const;
    double getSimTimePerSecond() const;

    void temperatureRescale(const double&);

    double getMeankT() const;
    double getMeanSqrkT() const;
    double getCurrentkT() const;

    Vector getMeanMomentum() const { return _sysMomentum.mean(); }
    Vector getCurrentMomentum() const { return _sysMomentum.current(); }

    double getTotalEnergy() const { return _internalE.current() + _KE.current(); }

    double getMeanUConfigurational() const;
    double getMeanSqrUConfigurational() const;
    inline double getConfigurationalU() const { return _internalE.current(); }

    Matrix getPressureTensor() const;

  protected:
    void stream(double dt);

    typedef std::pair<classKey, EEventType> CounterKey;

    struct CounterData
    {
      CounterData(): count(0), netimpulse({0,0,0}), netKEchange(0), netUchange(0) {}
      size_t count;
      Vector netimpulse;
      double netKEchange;
      double netUchange;
    };

    std::map<CounterKey, CounterData> _counters;
    std::chrono::system_clock::time_point _starttime;
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
    magnet::math::LogarithmicTimeCorrelator<double> _bulkVisc;
    magnet::math::LogarithmicTimeCorrelator<Vector> _crossVisc;
    std::vector<magnet::math::LogarithmicTimeCorrelator<Vector> > _thermalDiffusion;
    std::vector<magnet::math::LogarithmicTimeCorrelator<Vector> > _mutualDiffusion;
    std::vector<double> _internalEnergy;
    std::vector<double> _speciesMasses;
    std::vector<Vector> _speciesMomenta;
    double _systemMass;

    Matrix collisionalP;

    /*! \brief Flag to note when the late init is complete. */
    bool _lateInitComplete;
  };
}

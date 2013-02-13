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
#include <magnet/math/histogram.hpp>
#include <tr1/unordered_map>

namespace dynamo {
  class OPMisc;

  class OPIntEnergyHist: public OutputPlugin
  {
  public:
    OPIntEnergyHist(const dynamo::Simulation*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void eventUpdate(const IntEvent&, const PairEventData&);

    virtual void eventUpdate(const GlobalEvent&, const NEventData&);

    virtual void eventUpdate(const LocalEvent&, const NEventData&);

    virtual void eventUpdate(const System&, const NEventData&, const double&);

    virtual void output(magnet::xml::XmlStream&);

    virtual void changeSystem(OutputPlugin*);
  
    void operator<<(const magnet::xml::Node&);

    std::tr1::unordered_map<int, double> getImprovedW() const;
    inline double getBinWidth() const { return intEnergyHist.getBinWidth(); }
  protected:
    void stream(double);

    void ticker();

    magnet::math::HistogramWeighted<> intEnergyHist;
    shared_ptr<const OPMisc> _ptrOPMisc;
    double weight;
    double binwidth;

  };
}

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
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <complex>

namespace dynamo {
  class OPOrientationalOrder: public OPTicker
  {
  public:
    OPOrientationalOrder(const dynamo::SimData*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();
  
    virtual void output(magnet::xml::XmlStream&);

    virtual void operator<<(const magnet::xml::Node&);

  protected:
    typedef std::complex<double> ComplexNum;

    std::vector<ComplexNum> _history;
    Vector _axis;
    double _rg;
    size_t _nblistID;
  };
}

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
#include <magnet/math/histogram.hpp>
#include <vector>
#include <boost/algorithm/string.hpp>


namespace dynamo {
  class OPRadialDistribution: public OPTicker
  {
  static const size_t N_moments = 3;

  public:
    OPRadialDistribution(const dynamo::Simulation*, 
			 const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();
  
    virtual void output(magnet::xml::XmlStream&);

    void operator<<(const magnet::xml::Node&);

    std::vector<std::pair<double, double> > getgrdata(size_t species1ID, size_t species2ID) const;
    double getBinWidth() const { return binWidth; }
  protected:
    double binWidth;
    size_t length;
    unsigned long _sampleCount;
    double sample_energy; 
    double sample_energy_bin_width;
    
    std::vector<std::vector<std::vector<unsigned long> > > gr_accumulator; //The accumulator for the g(r)
    std::vector<double> moments; //The moments of the accumulator, a running sum of each.
  };
}

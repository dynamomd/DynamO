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
#include "fuzzy_array.hpp"

namespace xml
{
  class XmlStream;
}

class C1DHistogram
{
 public:
  C1DHistogram(double binwidth):
    data(binwidth),
    sampleCount(0)
    {};
  
  C1DHistogram():
    sampleCount(0)
    {};
  
  void addVal(const double& val)
    {
      ++data[val + 0.5 * data.binWidth];
      ++sampleCount;
    }
  
  typedef std::pair<const long, unsigned long> lv1pair;
  
  void outputHistogram(xml::XmlStream&, double) const;
  
  CFuzzyArray<unsigned long> data;
  
  unsigned long sampleCount;
};

class C1DWeightHistogram
{
 public:
  C1DWeightHistogram(double binwidth):
    data(binwidth),
    sampleCount(0.0)
    {};
  
  C1DWeightHistogram():
    sampleCount(0.0)
    {};
  
  void addVal(double val, double weight)
    {
      data[val + 0.5 * data.binWidth] += weight;
      sampleCount += weight;
    }
  
  void resetBinWidth(double val)
  {
    data = CFuzzyArray<double>(val);
  }

  typedef std::pair<const long, double> lv1pair;
  
  void outputHistogram(xml::XmlStream&, double) const;
  void outputClearHistogram(xml::XmlStream&, double) const;
  
  CFuzzyArray<double> data;
  
  double sampleCount;
};

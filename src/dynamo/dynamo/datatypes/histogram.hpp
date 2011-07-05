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
#include <magnet/containers/fuzzy_array.hpp>

namespace magnet { namespace xml { class XmlStream; } }

class C1DHistogram : public magnet::containers::FuzzyArray<unsigned long>
{
  typedef magnet::containers::FuzzyArray<unsigned long> Container;

public:
  C1DHistogram(double binwidth):
    Container(binwidth),
    sampleCount(0)
    {}
  
  C1DHistogram():
    sampleCount(0)
    {}
  
  void addVal(const double& val)
    {
      ++operator[](val);
      ++sampleCount;
    }
  
  void outputHistogram(magnet::xml::XmlStream&, double) const;

  inline unsigned long getSampleCount() const { return sampleCount; }

protected:
  unsigned long sampleCount;
};

class C1DWeightHistogram: public magnet::containers::FuzzyArray<double>
{
  typedef magnet::containers::FuzzyArray<double> Container;
public:
  C1DWeightHistogram(double binwidth):
    Container(binwidth),
    sampleCount(0)
    {}
  
  C1DWeightHistogram():
    sampleCount(0)
    {}
  
  void addVal(double val, double weight)
    {
      Container::operator[](val) += weight;
      sampleCount += weight;
    }
    
  void outputHistogram(magnet::xml::XmlStream&, double) const;
  void outputClearHistogram(magnet::xml::XmlStream&, double) const;
  
  inline double getSampleCount() const { return sampleCount; }
protected:
  double sampleCount;
};

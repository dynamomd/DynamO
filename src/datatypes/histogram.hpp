/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CHistogram_H
#define CHistogram_H

#include "fuzzy_array.hpp"

namespace xmlw
{
  class XmlStream;
}

class C1DHistogram
{
 public:
  C1DHistogram(Iflt binwidth):
    data(binwidth),
    sampleCount(0)
    {};
  
  C1DHistogram():
    sampleCount(0)
    {};
  
  void addVal(const Iflt& val)
    {
      ++data[val + 0.5 * data.binwidth];
      ++sampleCount;
    }
  
  typedef std::pair<const long, unsigned long> lv1pair;
  
  void outputHistogram(xmlw::XmlStream&, Iflt) const;
  
  CFuzzyArray<unsigned long> data;
  
  unsigned long sampleCount;
};

class C1DWeightHistogram
{
 public:
  C1DWeightHistogram(Iflt binwidth):
    data(binwidth),
    sampleCount(0.0)
    {};
  
  C1DWeightHistogram():
    sampleCount(0.0)
    {};
  
  void addVal(Iflt val, Iflt weight)
    {
      data[val + 0.5 * data.binwidth] += weight;
      sampleCount += weight;
    }
  
  void resetBinWidth(Iflt val)
  {
    data = CFuzzyArray<Iflt>(val);
  }

  typedef std::pair<const long, Iflt> lv1pair;
  
  void outputHistogram(xmlw::XmlStream&, Iflt) const;
  void outputClearHistogram(xmlw::XmlStream&, Iflt) const;
  
  CFuzzyArray<Iflt> data;
  
  Iflt sampleCount;
};

#endif


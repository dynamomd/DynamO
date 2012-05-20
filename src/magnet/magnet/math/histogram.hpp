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
#include <magnet/containers/fuzzy_array.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>

namespace magnet {
  namespace math {
    template<bool shiftBin=false>
    class Histogram : public magnet::containers::FuzzyArray<unsigned long, shiftBin>
    {
      typedef typename magnet::containers::FuzzyArray<unsigned long,shiftBin> Container;

    public:
      Histogram(double binwidth):
	Container(binwidth),
	sampleCount(0)
      {}
  
      Histogram():
	sampleCount(0)
      {}
  
      void addVal(const double& val)
      {
	++(Container::operator[](val));
	++sampleCount;
      }
  
      void outputHistogram(magnet::xml::XmlStream& XML, double scalex) const
      {

	XML << magnet::xml::tag("Histogram")
	    << magnet::xml::attr("SampleCount")
	    << sampleCount
	    << magnet::xml::attr("Dimension") << 1
	    << magnet::xml::attr("BinWidth") << Container::getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const typename Container::value_type &p1, *this)
	  avgSum += (static_cast<double>(p1.first) + 0.5 * shiftBin) * p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << (avgSum * Container::getBinWidth() * scalex / sampleCount)
	    << magnet::xml::chardata();
  
	BOOST_FOREACH(const typename Container::value_type &p1, *this)
	  XML << (p1.first + 0.5 * shiftBin) * Container::getBinWidth() * scalex << " " 
	      << static_cast<double>(p1.second)
	  /(Container::getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("Histogram");
      }

      inline unsigned long getSampleCount() const { return sampleCount; }

    protected:
      unsigned long sampleCount;
    };

    template<bool shiftBin=false>
    class HistogramWeighted: public magnet::containers::FuzzyArray<double,shiftBin>
    {
      typedef magnet::containers::FuzzyArray<double> Container;
    public:
      HistogramWeighted(double binwidth):
	Container(binwidth),
	sampleCount(0)
      {}
  
      HistogramWeighted():
	sampleCount(0)
      {}
  
      void addVal(double val, double weight)
      {
	Container::operator[](val) += weight;
	sampleCount += weight;
      }
    
      void outputHistogram(magnet::xml::XmlStream& XML, double scalex) const
      {
	XML << magnet::xml::tag("HistogramWeighted")
	    << magnet::xml::attr("TotalWeight")
	    << sampleCount
	    << magnet::xml::attr("Dimension") << 1
	    << magnet::xml::attr("BinWidth") << Container::getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  avgSum += static_cast<double>(p1.first + 0.5 * shiftBin) * p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << avgSum * Container::getBinWidth() * scalex / sampleCount
	    << magnet::xml::chardata();
  
	//This gives mathmatically correct but not really pretty
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  XML << (p1.first + 0.5 * shiftBin) * Container::getBinWidth() * scalex << " "
	      << static_cast<double>(p1.second)
	  / (Container::getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("HistogramWeighted");
      }


      void outputClearHistogram(magnet::xml::XmlStream& XML, 
				double scalex) const
      {
	XML << magnet::xml::tag("HistogramWeighted")
	    << magnet::xml::attr("TotalWeight")
	    << sampleCount
	    << magnet::xml::attr("Dimension") << 1
	    << magnet::xml::attr("BinWidth") << Container::getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  avgSum += static_cast<double>(p1.first + 0.5 * shiftBin)* p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << (avgSum * Container::getBinWidth() / sampleCount)
	    << magnet::xml::chardata();
    
	//This one gives histograms usable by the reweight program
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  XML << (p1.first + 0.5 * shiftBin) * Container::getBinWidth() * scalex << " " 
	      << static_cast<double>(p1.second)
	  / (Container::getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("HistogramWeighted");
      }
  
      inline double getSampleCount() const { return sampleCount; }
    protected:
      double sampleCount;
    };
  }
}

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
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>

namespace magnet {
  namespace math {
    class Histogram : public magnet::containers::FuzzyArray<unsigned long>
    {
      typedef magnet::containers::FuzzyArray<unsigned long> Container;

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
	++operator[](val);
	++sampleCount;
      }
  
      void outputHistogram(magnet::xml::XmlStream& XML, double scalex) const
      {

	XML << magnet::xml::tag("Histogram")
	    << magnet::xml::attr("SampleCount")
	    << sampleCount
	    << magnet::xml::attr("Dimension") << 1
	    << magnet::xml::attr("BinWidth") << getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  avgSum += (static_cast<double>(p1.first) + 0.5) * p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << (avgSum * getBinWidth() * scalex / sampleCount)
	    << magnet::xml::chardata();
  
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  XML << p1.first * getBinWidth() * scalex << " " 
	      << static_cast<double>(p1.second)
	  /(getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("Histogram");
      }

      inline unsigned long getSampleCount() const { return sampleCount; }

    protected:
      unsigned long sampleCount;
    };

    class HistogramWeighted: public magnet::containers::FuzzyArray<double>
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
	    << magnet::xml::attr("BinWidth") << getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  avgSum += static_cast<double>(p1.first) * p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << avgSum * getBinWidth() * scalex / sampleCount
	    << magnet::xml::chardata();
  
	//This gives mathmatically correct but not really pretty
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  XML << p1.first * getBinWidth() * scalex << " "
	      << static_cast<double>(p1.second)
	  / (getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("HistogramWeighted");
      }


      void outputClearHistogram(magnet::xml::XmlStream& XML, 
				double scalex) const
      {
	XML << magnet::xml::tag("HistogramWeighted")
	    << magnet::xml::attr("TotalWeight")
	    << sampleCount
	    << magnet::xml::attr("Dimension") << 1
	    << magnet::xml::attr("BinWidth") << getBinWidth() * scalex;
  
	double avgSum = 0.0;
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  avgSum += static_cast<double>(p1.first)* p1.second;
  
	XML << magnet::xml::attr("AverageVal")
	    << (avgSum * getBinWidth() / sampleCount)
	    << magnet::xml::chardata();
    
	//This one gives histograms usable by the reweight program
	BOOST_FOREACH(const Container::value_type &p1, *this)
	  XML << p1.first * getBinWidth() * scalex << " " 
	      << static_cast<double>(p1.second)
	  / (getBinWidth() * sampleCount * scalex) << "\n";
  
	XML << magnet::xml::endtag("HistogramWeighted");
      }
  
      inline double getSampleCount() const { return sampleCount; }
    protected:
      double sampleCount;
    };
  }
}

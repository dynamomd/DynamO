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

#include "histogram.hpp"
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

void 
C1DHistogram::outputHistogram(magnet::xml::XmlStream& XML, double scalex) const
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
  
///////Pretty histogram output method
//  long lastx = data.data.begin()->first - 1;
//
//  XML << data.data.begin()->first * data.binWidth * scalex 
//      << " " << 0 << "\n";
//
//      BOOST_FOREACH(const lv1pair &p1, data.data)
//	{
//	  double y = static_cast<double>(p1.second)
//	    /(data.binWidth * sampleCount * scalex);
//	  
//	  if (p1.first - 1 != lastx)
//	    XML << (lastx + 1) * data.binWidth * scalex  << " " << 0 << "\n"
//		<< p1.first * data.binWidth * scalex << " " << 0 << "\n";
//	  
//	  XML << p1.first * data.binWidth * scalex << " " 
//	      << y << "\n"
//	      << (p1.first + 1) * data.binWidth * scalex 
//	      << " " << y << "\n";
//	  
//	  lastx = p1.first;
//	}
//      
//      XML << (lastx + 1) * data.binWidth * scalex
//	  << " " << 0 << "\n";

  BOOST_FOREACH(const Container::value_type &p1, *this)
    XML << p1.first * getBinWidth() * scalex << " " 
	<< static_cast<double>(p1.second)
    /(getBinWidth() * sampleCount * scalex) << "\n";
  
  XML << magnet::xml::endtag("Histogram");
}

void 
C1DWeightHistogram::outputHistogram(magnet::xml::XmlStream & XML, double scalex) const
{
  XML << magnet::xml::tag("WeightHistogram")
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
  

  //This gives pretty but not really useful drawings of histograms
      //long lastx = data.data.begin()->first - 1;
      //
      //XML << data.data.begin()->first * data.binWidth * scalex 
      //	  << " " << 0 << "\n";
      //
      //BOOST_FOREACH(const lv1pair &p1, data.data)
      //	{
      //	  double y = static_cast<double>(p1.second)
      //	    /(data.binWidth * sampleCount * scalex);
      //	  
      //	  if (p1.first - 1 != lastx)
      //	    XML << (lastx + 1) * data.binWidth * scalex  << " " << 0 << "\n"
      //		<< p1.first * data.binWidth * scalex << " " << 0 << "\n";
      //	  
      //	  XML << p1.first * data.binWidth * scalex << " " 
      //	      << y << "\n"
      //	      << (p1.first + 1) * data.binWidth * scalex 
      //	      << " " << y << "\n";
      //	  
      //	  lastx = p1.first;
      //	}
      //
      //XML << lastx * data.binWidth * scalex
      //	  << " " << 0 << "\n";
      //

      //This gives mathmatically correct but not really pretty
  BOOST_FOREACH(const Container::value_type &p1, *this)
    XML << p1.first * getBinWidth() * scalex << " "
	<< static_cast<double>(p1.second)
    / (getBinWidth() * sampleCount * scalex) << "\n";
  
  XML << magnet::xml::endtag("WeightHistogram");
}

void 
C1DWeightHistogram::outputClearHistogram(magnet::xml::XmlStream & XML, double scalex) const
{
  XML << magnet::xml::tag("WeightHistogram")
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
  
  XML << magnet::xml::endtag("WeightHistogram");
}

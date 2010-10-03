/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../extcode/xmlwriter.hpp"

void 
C1DHistogram::outputHistogram(xml::XmlStream& XML, double scalex) const
{

  XML << xml::tag("Histogram")
      << xml::attr("SampleCount")
      << sampleCount
      << xml::attr("Dimension") << 1
      << xml::attr("BinWidth") << data.binWidth * scalex;
  
  double avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += (static_cast<double>(p1.first) + 0.5) * p1.second;
  
  XML << xml::attr("AverageVal")
      << (avgSum * data.binWidth * scalex / sampleCount)
      << xml::chardata();
  
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

  BOOST_FOREACH(const lv1pair &p1, data.data)
    XML << p1.first * data.binWidth * scalex << " " 
	<< static_cast<double>(p1.second)
    /(data.binWidth * sampleCount * scalex) << "\n";
  
  XML << xml::endtag("Histogram");
}

void 
C1DWeightHistogram::outputHistogram(xml::XmlStream & XML, double scalex) const
{
  XML << xml::tag("WeightHistogram")
      << xml::attr("TotalWeight")
      << sampleCount
      << xml::attr("Dimension") << 1
      << xml::attr("BinWidth") << data.binWidth * scalex;
  
  double avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += static_cast<double>(p1.first) * p1.second;
  
  XML << xml::attr("AverageVal")
      << avgSum * data.binWidth * scalex / sampleCount
      << xml::chardata();
  

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
  BOOST_FOREACH(const lv1pair &p1, data.data)
    XML << p1.first * data.binWidth * scalex << " "
	<< static_cast<double>(p1.second)
    / (data.binWidth * sampleCount * scalex) << "\n";
  
  XML << xml::endtag("WeightHistogram");
}

void 
C1DWeightHistogram::outputClearHistogram(xml::XmlStream & XML, double scalex) const
{
  XML << xml::tag("WeightHistogram")
      << xml::attr("TotalWeight")
      << sampleCount
      << xml::attr("Dimension") << 1
      << xml::attr("BinWidth") << data.binWidth * scalex;
  
  double avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += static_cast<double>(p1.first)* p1.second;
  
  XML << xml::attr("AverageVal")
      << (avgSum * data.binWidth / sampleCount)
      << xml::chardata();
    
  //This one gives histograms usable by the reweight program
  BOOST_FOREACH(const lv1pair &p1, data.data)
    XML << p1.first * data.binWidth * scalex << " " 
	<< static_cast<double>(p1.second)
    / (data.binWidth * sampleCount * scalex) << "\n";
  
  XML << xml::endtag("WeightHistogram");
}

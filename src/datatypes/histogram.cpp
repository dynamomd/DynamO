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

#include "histogram.hpp"
#include <boost/foreach.hpp>
#include "../extcode/xmlwriter.hpp"

void 
C1DHistogram::outputHistogram(xmlw::XmlStream& XML, Iflt scalex) const
{

  XML << xmlw::tag("Histogram")
      << xmlw::attr("SampleCount")
      << sampleCount
      << xmlw::attr("Dimension") << 1;
  
  Iflt avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += (static_cast<Iflt>(p1.first) + 0.5) * p1.second;
  
  XML << xmlw::attr("AverageVal")
      << (avgSum * data.binWidth * scalex / sampleCount)
      << xmlw::chardata();
  
  
  long lastx = data.data.begin()->first - 1;

  XML << data.data.begin()->first * data.binWidth * scalex 
      << " " << 0 << "\n";

  BOOST_FOREACH(const lv1pair &p1, data.data)
    {
      Iflt y = static_cast<Iflt>(p1.second)
	/(data.binWidth * sampleCount * scalex);

      if (p1.first - 1 != lastx)
	XML << (lastx + 1) * data.binWidth * scalex  << " " << 0 << "\n"
	    << p1.first * data.binWidth * scalex << " " << 0 << "\n";
      
      XML << p1.first * data.binWidth * scalex << " " 
	  << y << "\n"
	  << (p1.first + 1) * data.binWidth * scalex 
	  << " " << y << "\n";

      lastx = p1.first;
    }
  
  XML << (lastx + 1) * data.binWidth * scalex
      << " " << 0 << "\n";
  
  XML << xmlw::endtag("Histogram");
}

void 
C1DWeightHistogram::outputHistogram(xmlw::XmlStream & XML, Iflt scalex) const
{
  XML << xmlw::tag("WeightHistogram")
      << xmlw::attr("TotalWeight")
      << sampleCount
      << xmlw::attr("Dimension") << 1;
  
  Iflt avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += (static_cast<Iflt>(p1.first) + 0.5) * p1.second;
  
  XML << xmlw::attr("AverageVal")
      << (avgSum * data.binWidth * scalex / sampleCount)
      << xmlw::chardata();
  

  if (false)
    {
      //This gives pretty but not really useful drawings of histograms
      long lastx = data.data.begin()->first - 1;
      
      XML << data.data.begin()->first * data.binWidth * scalex 
	  << " " << 0 << "\n";
      
      BOOST_FOREACH(const lv1pair &p1, data.data)
	{
	  Iflt y = static_cast<Iflt>(p1.second)
	    /(data.binWidth * sampleCount * scalex);
	  
	  if (p1.first - 1 != lastx)
	    XML << (lastx + 1) * data.binWidth * scalex  << " " << 0 << "\n"
		<< p1.first * data.binWidth * scalex << " " << 0 << "\n";
	  
	  XML << p1.first * data.binWidth * scalex << " " 
	      << y << "\n"
	      << (p1.first + 1) * data.binWidth * scalex 
	      << " " << y << "\n";
	  
	  lastx = p1.first;
	}

      XML << lastx * data.binWidth * scalex
	  << " " << 0 << "\n";

    } else
    {
      //This gives mathmatically correct but not really pretty
      BOOST_FOREACH(const lv1pair &p1, data.data)
	XML << (p1.first + 0.5) * data.binWidth * scalex << " "
	    << static_cast<Iflt>(p1.second)
	/ (data.binWidth * sampleCount * scalex) << "\n";
    }

  
  XML << xmlw::endtag("WeightHistogram");
}

void 
C1DWeightHistogram::outputClearHistogram(xmlw::XmlStream & XML, Iflt scalex) const
{
  XML << xmlw::tag("WeightHistogram")
      << xmlw::attr("TotalWeight")
      << sampleCount
      << xmlw::attr("Dimension") << 1;
  
  Iflt avgSum = 0.0;
  BOOST_FOREACH(const lv1pair &p1, data.data)
    avgSum += (static_cast<Iflt>(p1.first) + 0.5) * p1.second;
  
  XML << xmlw::attr("AverageVal")
      << (avgSum * data.binWidth / sampleCount)
      << xmlw::chardata();
    
  //This one gives histograms usable by the reweight program
  BOOST_FOREACH(const lv1pair &p1, data.data)
    XML << (p1.first + 0.5) * data.binWidth * scalex << " " 
	<< static_cast<Iflt>(p1.second)
    / (data.binWidth * sampleCount * scalex) << "\n";
  
  XML << xmlw::endtag("WeightHistogram");
}

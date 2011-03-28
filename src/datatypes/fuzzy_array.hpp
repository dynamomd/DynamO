/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include "../base/constants.hpp"
#include <magnet/exception.hpp>
#include <vector>
#include <map>
#include <cmath>

template<class T>
class CFuzzyArray
{
public:  
  typedef std::pair<const long, T> mapType;

  CFuzzyArray(double binwidth):
    binWidth(binwidth)
  {}
  
  CFuzzyArray():
    binWidth(0.0)
  {}

  void setBinWidth(double bw)
    {
      binWidth = bw;
      data.clear();
    }
  
  T& operator[](const double& x)
  { 
    double val = x/binWidth;
    return data[static_cast<long>(val) - static_cast<long>(val < 0)]; 
  }

  inline T &operator[](const long &x)
  {
    return data[x];
  }
  
  double binWidth;
  
  std::map<long, T> data;
};

template<class T>
class CFuzzyArray<CFuzzyArray<T> >
{
public:
  typedef std::pair<const long, CFuzzyArray<T> > mapType;
  
  CFuzzyArray(double binwidth):
    binWidth(binwidth)
  {}
  
  CFuzzyArray():
    binWidth(0.0)
  {}

  void setBinWidth(double bw)
    {
      binWidth = bw;
      data.clear();
    }
  

  CFuzzyArray<T>& operator[](const double&x)
    {
      long i = lround((x/binWidth)-0.5);
      
      if (data.find(i) == data.end())
	data[i].setBinWidth(binWidth);
      
      return data[i];
    }

  CFuzzyArray<T> &operator[](const long &i)
  {
    if (data.find(i) == data.end())
      data[i].setBinWidth(binWidth);

    return data[i];
  }
  
  double binWidth;
  
  std::map<long, CFuzzyArray<T> > data;
};


template<class T>
class CFuzzyArray2
{
public:  
  CFuzzyArray2(double binwidth, double origin2, long nbins):
    binWidth(binwidth),
    origin(origin2),
    data(nbins, T(0))
    {}
  
    T& operator[](const double &x)
      {
	long i = static_cast<long>((x-origin)/binWidth);
	if (i > static_cast<long>(data.size()))
	  M_throw() << "Data too high, " << i;
	if (i < 0)
	  M_throw() << "Data too low, " << i;
	
	return data[i];
      }
    
    T& operator[](const long &x)
      {
	if (x > static_cast<long>(data.size()))
	  M_throw() << "Data too high, " << x;
	if (x < 0)
	  M_throw() << "Data too low, " << x;
	
	return data[x];
      }
  
  double binWidth, origin;
  std::vector<T> data;
};

template<class T>
class CFuzzyArray2<CFuzzyArray2<T> >
{
 public:
  
 CFuzzyArray2(double binwidth, double origin2, long nbins):
  binWidth(binwidth),
    origin(origin2),
    data(nbins, CFuzzyArray2<T>(binwidth,origin2,nbins))
      {}
  
  CFuzzyArray2<T>& operator[](const double &x)
    {
      long i = static_cast<long>((x-origin)/binWidth);
      if (i > static_cast<long>(data.size()))
	M_throw() << "Data too high, " << i;
      if (i < 0)
	M_throw() << "Data too low, " << i;
      return data[i];
    }
  
  CFuzzyArray2<T>& operator[](const long &x)
    {
      if (x > static_cast<long>(data.size()))
	M_throw() << "Data too high, " << x;
      if (x < 0)
	M_throw() << "Data too low, " << x;
      
      return data[x];
    }
  
  double binWidth, origin;
  std::vector<CFuzzyArray2<T> > data;
};

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
#include <magnet/exception.hpp>
#include <map>
#include <cmath>

namespace magnet {
  namespace containers {
    /*! \brief A class to create a sparse array which is addressed
      using floating point values.
     
      Standard arrays/vectors are addressed using a unsigned integer
      type, specifying the offset into the data. However, there are
      times where it would be easier to address an array using a
      floating point type (the prime example is in
      histogramming). This class solves this problem by storing a bin
      width (actually the inverse bin width \ref _invBinWidth), which
      is used to map a floating point value to a bin.
     
      This class is space efficient as it uses a map to store the
      allocated bins. A map (not a unordered_map) must be the
      default, as histogramming output assumes the values are sorted.
     
      Both the bin width and the inverse bin width are stored within
      the class. The inverse is stored as multiply operations are far
      cheaper than divides, so this accelerates the array access
      operator. The bin width is stored to ensure the following code
      is true
      \code
      FuzzyArray<T>(some_bin_width).getBinWidth() == some_bin_width
      \endcode

      There may be a loss of precision when converting the inverse bin
      width back to the bin width which would cause the above
      comparison to fail. As such, the original bin width is stored.

      \tparam T The type stored by the FuzzyArray.
      \tparam shiftBin If false, bins are centered about whole integers, if true bins are centered between whole integers.
      \tparam Container The underlying container used in the FuzzyArray.
    */
    template<class T, bool shiftBin = false, class Container = std::map<long, T> >
    class FuzzyArray : public Container
    {
    public:
      /*! \brief Default constructor.
        
	\param binwidth The width of a single bin.
       */
      FuzzyArray(double binwidth = 1):
	_binWidth(binwidth),
	_invBinWidth(1 / binwidth)
      {}
      
      /*! \brief Set the width of the bins of the FuzzyArray.
       
        \note This clears all stored values in the FuzzyArray.
       */
      void setBinWidth(double bw)
      {
	_binWidth = bw;
	_invBinWidth = 1 / bw;
	Container::clear();
      }
     
      /*! \brief Returns the current bin width. */
      double getBinWidth() const
      { return _binWidth; }
  
      /*! \brief Access an element of the FuzzyArray.
       
        If an addressed element does not exist in the underlying
        container, it is initialized using the default constructor
        T(). This initializes POD types (e.g., float, double, int) to
        0.
       */
      T& operator[](const double& x)
      { return Container::operator[](lrint(x * _invBinWidth + shiftBin * 0.5)); }

    protected:
      double _binWidth;
      double _invBinWidth;
    };
  }
}

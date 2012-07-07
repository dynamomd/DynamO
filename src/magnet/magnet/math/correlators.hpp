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
#include <boost/circular_buffer.hpp>
#include <vector>
#include <utility>

namespace magnet {
  namespace math {
    /*! \brief A class for taking the Einstein correlation of a series of
      values.

      We often want to collect functions in the following form
      \f[f_{1,2}(j)=\left\langle \left[W^{(1)}_{i} \cdot
      W^{(1)}_{i+j}\right]\right\rangle_i\f]
    
      Such functions are Einstein correlators and the values
      \f$W^{(1)}_{i}\f$ and \f$W^{(2)}_{i}\f$ are the integrals of the
      microscopic fluxes. The angle brackets \f$\langle\rangle_i\f$
      indicate an average over the origins \f$i\f$.

      \tparam T The type of variable to be correlated. This must have a
      default constructor which zeros the type.
    */
    template<class T, size_t length>
    class Correlator
    {
    public:
      Correlator() { clear(); }

      /*! \brief Add a new pair of
       */
      void push(const T& val1, const T& val2)
      { 
	_sample_history.push_front(std::pair<T, T>(val1, val2));
	pass(); 
      }

      void clear()
      {
	_count = 0;

	_sample_history.clear();
	_sample_history.set_capacity(length);

	_correlator.clear();
	//In the next line, we rely on the default constructed type being
	//whatever passes for zero. This is true for 
	_correlator.resize(length);
      }

      std::vector<T> getAveragedCorrelator()
      {
	std::vector<T> avg_correlator;
	avg_correlator.resize(_sample_history.size());
    
	for (size_t i(0); i < avg_correlator.size(); ++i)
	  avg_correlator[i] = _correlator[i] / (_count - i);
    
	return avg_correlator;
      }

    protected:
      boost::circular_buffer<std::pair<T, T> > _sample_history;
      std::vector<T> _correlator;
      size_t _count;

      void pass()
      {
	++_count;
	size_t sample_length = std::min(length, _count);

	//Here, we're using the default constructor again to make a sum
	//variable with a value of zero.
	std::pair<T,T> sum = std::pair<T,T>();
	for (size_t i(0); i < _sample_history.size(); ++i)
	  {
	    sum.first += _sample_history[i].first;
	    sum.second += _sample_history[i].second;
	    _correlator[i] += sum.first * sum.second;
	  }
      }
    };

    template<class T, size_t length>
    class TimeCorrelator: protected Correlator<T,length>
    {
      typedef Correlator<T,length> Base;
    public:
      TimeCorrelator(double sample_time):_sample_time(sample_time) { clear(); }

      void impulse(const T& val1, const T& val2)
      {
	_impulse_values.first += val1; 
	_impulse_values.second += val2;
      }

      void freestream(double dt)
      {

      }

      void clear()
      {
	Base::clear();
	_freestream_values = std::pair<T,T>();
	_impulse_values = std::pair<T,T>();
	_sample_time = _current_time;
      }

      using Base::getAveragedCorrelator;

    protected:
      std::pair<T,T> _freestream_values;
      std::pair<T,T> _impulse_values;
      double _sample_time;
      double _current_time;
    };
  }
}


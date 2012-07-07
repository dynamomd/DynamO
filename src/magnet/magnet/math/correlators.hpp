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
#include <boost/foreach.hpp>
#include <vector>
#include <utility>
#include <exception>
#include <tr1/tuple>

namespace magnet {
  namespace math {
    /*! \brief A class for taking the Einstein correlation of a series of
      values.

      We often want to collect functions in the following form
      \f[f_{1,2}(j)=\left\langle \left[W^{(1)}_{i} \cdot
      W^{(2)}_{i+j}\right]\right\rangle_i\f]
    
      Such functions are Einstein correlators and the values
      \f$W^{(1)}_{i}\f$ and \f$W^{(2)}_{i}\f$ are the integrals of the
      microscopic fluxes over some interval of length \f$j\f$ and
      origin \f$i\f$. The angle brackets \f$\langle\rangle_i\f$
      indicate an average over the origins \f$i\f$.

      The correlators are calculated as described in "Molecular
      Dynamics Simulation: Elementary Methods," by J. M. Haile. We
      store a list of differences of the \f$W\f$ values between each
      new origin (\f$\Delta W^{(1)}_i=W^{(1)}_{i+1} -
      W^{(1)}_{i}\f$). This allows us to have a rolling correlation
      window and to collect the maximum amount of data. 

      Typically, we have to integrate the free streaming and impulsive
      contributions to the differences of \f$W\f$, and this
      functionality is provided by the \ref TimeCorrelator class.

      \tparam T The type of variable to be correlated, allowing vector
      and matrix values to be correlated. This must have a default
      constructor which zeros the type.
    */
    template<class T>
    class Correlator
    {
    public:
      /* \brief Create a Correlator with a given length (i.e., maximum
	 value of \f$j\f$ calculable for \f$f_{1,2}(j)\f$
	 
	 \param length The maximum length of the correlator.
       */
      Correlator(size_t length): _length(length) { clear(); }

      /*! \brief Add a new pair of \f$\Delta W^{(1)}\f$ and \f$\Delta
          W^{(2)}\f$ values to the correlator

	  \param W1 \f$\Delta W^{(1)}\f$
	  \param W2 \f$\Delta W^{(2)}\f$
       */
      void push(const T& W1, const T& W2)
      { 
	_sample_history.push_front(std::pair<T, T>(W1, W2));
	pass(); 
      }

      /*! \brief Clear the correlator. */
      void clear()
      {
	_count = 0;

	_sample_history.clear();
	_sample_history.set_capacity(_length);

	_correlator.clear();
	//In the next line, we rely on the default constructed type
	//being whatever passes for zero. This is true for all
	//built-in types.
	_correlator.resize(_length);
      }
      
      /*! \brief Returns a vector where each component is a different
	value of \f$j\f$ for the function \f$f_{1,2}(j)=\left\langle
	\left[W^{(1)}_{i} \cdot W^{(2)}_{i+j}\right]\right\rangle_i\f$.
	
	This list of values may be less than the length of the
	correlator if less \f$\Delta W\f$ have been push()'ed than the
	length of the correlator.
      */
      std::vector<T> getAveragedCorrelator()
      {
	std::vector<T> avg_correlator;
	avg_correlator.resize(_sample_history.size());
    
	for (size_t i(0); i < avg_correlator.size(); ++i)
	  avg_correlator[i] = _correlator[i] / (_count - i);
    
	return avg_correlator;
      }

      /*! \brief Returns the number of samples collected for the
          \f$j\f$th correlation of \f$f_{1,2}(j)\f$. 
      */
      size_t getSampleCount(size_t i) const
      {
	return _count - i;
      }

    protected:
      boost::circular_buffer<std::pair<T, T> > _sample_history;
      std::vector<T> _correlator;
      size_t _count;
      size_t _length;

      /*! \brief Performs a pass of the correlator, calculating all of
          the \f$\left[W^{(1)}_{i} \cdot W^{(2)}_{i+j}\right]\f$
          values it can with the current data. 

	  These are summed up, ready to be divided by the pass count
	  in getAveragedCorrelator() to return averaged values.
      */
      void pass()
      {
	++_count;
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

    /*! \brief A modification of the Correlator class to perform the
        integration of the
     */
    template<class T>
    class TimeCorrelator: protected Correlator<T>
    {
      typedef Correlator<T> Base;
    public:
      TimeCorrelator(double sample_time, size_t length): 
	Base(length),
	_sample_time(sample_time)
      { 
	if (sample_time <= 0)
	  throw std::runtime_error("TimeCorrelator requires a positive, non-zero sample time");

	clear(); 
      }

      void addImpulse(const T& val1, const T& val2)
      {
	_W_sums.first += val1; 
	_W_sums.second += val2;
      }

      void setFreeStreamValue(const T& val1, const T& val2)
      {
	_freestream_values = std::pair<T,T>(val1,val2);
      }

      void freeStream(double dt)
      {
	while ((_current_time + dt) >= _sample_time)
	  {
	    const double deltat = _sample_time - _current_time;
	    _W_sums.first += _freestream_values.first * deltat;
	    _W_sums.second += _freestream_values.second * deltat;
	    Base::push(_W_sums.first, _W_sums.second);

	    _W_sums.first = _W_sums.second = _current_time = 0;
	    dt -= deltat;
	  }

	_W_sums.first += _freestream_values.first * dt;
	_W_sums.second += _freestream_values.second * dt;
	_current_time += dt;
      }

      void clear()
      {
	Base::clear();
	_freestream_values = std::pair<T,T>();
	_W_sums = std::pair<T,T>();
	_current_time = 0;
	
      }

      using Base::getAveragedCorrelator;
      using Base::getSampleCount;

      double getSampleTime() const { return _sample_time; }

    protected:
      std::pair<T,T> _freestream_values;
      std::pair<T,T> _W_sums;
      double _sample_time;
      double _current_time;
    };

    template<class T>
    class LogarithmicTimeCorrelator
    {
      typedef std::vector<TimeCorrelator<T> > Container;
    public:
      LogarithmicTimeCorrelator(double sample_time, size_t length, size_t scaling = 2):
	_sample_time(sample_time),
	_length(length),
	_scaling(scaling)
      {
	clear();
      }
    
      void clear()
      {
	_current_time = 0;
	_freestream_values = _freestream_sum = _impulse_sum= std::pair<T,T>();
	_sample_time /= (1 << _correlators.size());
	_correlators.clear();
      }

      void addImpulse(const T& val1, const T& val2)
      {
	_impulse_sum.first += val1; 
	_impulse_sum.second += val2;

	BOOST_FOREACH(TimeCorrelator<T>& correlator, _correlators)
	  correlator.addImpulse(val1, val2);
      }

      void setFreeStreamValue(const T& val1, const T& val2)
      {
	_freestream_values = std::pair<T,T>(val1, val2);

	BOOST_FOREACH(TimeCorrelator<T>& correlator, _correlators)
	  correlator.setFreeStreamValue(val1, val2);
      }

      void freeStream(const double dt)
      {
	//Check if we need to add a new correlator
	while ((_current_time + dt) >= _sample_time)
	  {
	    //Add a new correlator
	    _correlators.push_back(TimeCorrelator<T>(_sample_time, _length));
	    _sample_time *= _scaling;
	    
	    //Pretend the correlator has actually been here all along, gathering impulse data
	    TimeCorrelator<T>& new_correlator = _correlators.back();
	    new_correlator.addImpulse(_impulse_sum.first, _impulse_sum.second);

	    //Also fake a freestream integration, it is fine to fake
	    //this as the correlators cant resolve smaller than a
	    //_sample_time.
	    new_correlator.setFreeStreamValue(_freestream_sum.first / _current_time, _freestream_sum.second / _current_time);
	    new_correlator.freeStream(_current_time);

	    //Set the correct free streaming value
	    new_correlator.setFreeStreamValue(_freestream_values.first, _freestream_values.second);
	  }

	BOOST_FOREACH(TimeCorrelator<T>& correlator, _correlators)
	  correlator.freeStream(dt);

	_freestream_sum.first += _freestream_values.first * dt;
	_freestream_sum.second += _freestream_values.second * dt;
	_current_time += dt;
      }

      struct Data
      {
	Data(double t, size_t sc, T v): time(t), sample_count(sc), value(v) {}

	double time;
	size_t sample_count;
	T value;
      };

      std::vector<Data> getAveragedCorrelator()
      {
	std::vector<Data> avg_correlator;

	if (!_correlators.empty())
	  {
	    {
	      std::vector<T> result = _correlators.front().getAveragedCorrelator();
	      
	      for (size_t i(0); i < result.size(); ++i)
		avg_correlator.push_back(Data(_correlators.front().getSampleTime() * (i+1), 
					      _correlators.front().getSampleCount(i), result[i]));
	    }

	    //Now copy the rest of the correlators
	    for (size_t i(1); i < _correlators.size(); ++i)
	      {
		std::vector<T> result = _correlators[i].getAveragedCorrelator();
		for (size_t j(_length / _scaling); j < result.size(); ++j) 
		  avg_correlator.push_back(Data(_correlators[i].getSampleTime() * (j+1),
						_correlators[i].getSampleCount(j), result[j]));
	      }
	  }
	return avg_correlator;
      }


    protected:
      double _sample_time;
      double _current_time;
      size_t _length;
      size_t _scaling;
      std::pair<T,T> _freestream_values;
      std::pair<T,T> _impulse_sum;
      std::pair<T,T> _freestream_sum;
      
      std::vector<TimeCorrelator<T> > _correlators;
    };
  }
}


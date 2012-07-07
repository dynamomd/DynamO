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
#include <magnet/math/vector.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <utility>
#include <exception>
#include <tr1/tuple>

namespace magnet {
  namespace math {
    /*! \brief This class facilitates the exact averaging of
        properties during an event driven simulation.

	As many simulation properties are constant between events, we
	can actually calculate exact statistics on the fluctuations of
	these properties. This class calculates the mean, square mean,
	minimum and maximum of a given property.
     */
    class TimeAveragedProperty
    {
    public:
      TimeAveragedProperty():
	_current_value(0),
	_zero_moment(0),
	_first_moment(0),
	_second_moment(0),
	_min(HUGE_VAL),
	_max(-HUGE_VAL)
      {}

      void swapCurrentValues(TimeAveragedProperty& op)
      { std::swap(_current_value, op._current_value); }
	
      TimeAveragedProperty& operator=(double value)
      { 
	_current_value = value; 
	_min = std::min(_current_value, _min);
	_max = std::max(_current_value, _max);
	return *this;
      }

      TimeAveragedProperty& operator+=(double change)
      { return operator=(_current_value + change); }

      void stream(double dt)
      {
	_zero_moment += dt;
	_first_moment += dt * _current_value;
	_second_moment += dt * _current_value * _current_value;
      }

      double time() const { return _zero_moment; }

      double mean() const { return _first_moment / _zero_moment; }

      double meanSqr() const { return _second_moment / _zero_moment; }

      double min() const { return _min; }

      double max() const { return _max; }
      
      double current() const { return _current_value; }

      bool empty() const { return _zero_moment; }

    protected:
      double _current_value;
      double _zero_moment;
      double _first_moment;
      double _second_moment;
      double _min;
      double _max;
    };
  }
}


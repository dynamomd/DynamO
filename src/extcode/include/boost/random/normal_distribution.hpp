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

/* boost random/normal_distribution.hpp header file
 *
 * Copyright Jens Maurer 2000-2001
 * Copyright Marcus Bannerman 2007
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 * $Id: normal_distribution.hpp,v 1.20 2004/07/27 03:43:32 dgregor Exp $
 *
 * Revision history
 *  2001-02-18  moved to individual header files
 *  2007-11-10  M.B. Changed to the outward cartesian form of the Box-Muller 
 *              transform. 
 */

#ifndef BOOST_RANDOM_NORMAL_DISTRIBUTION_HPP
#define BOOST_RANDOM_NORMAL_DISTRIBUTION_HPP

#include <cmath>
#include <cassert>
#include <iostream>
#include <boost/limits.hpp>
#include <boost/static_assert.hpp>

namespace boost {

// deterministic polar method, uses trigonometric functions
template<class RealType = double>
class normal_distribution
{
public:
  typedef RealType input_type;
  typedef RealType result_type;

#if !defined(BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS) && !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300)
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif

  explicit normal_distribution(const result_type& mean = result_type(0),
                               const result_type& sigma = result_type(1))
    : _mean(mean), _sigma(sigma), _valid(false)
  {
    assert(sigma >= result_type(0));
  }

  // compiler-generated copy constructor is NOT fine, need to purge cache
  normal_distribution(const normal_distribution& other)
    : _mean(other._mean), _sigma(other._sigma), _r1(NAN), _r2(NAN), _valid(false)
  {}

  // compiler-generated copy ctor and assignment operator are fine

  RealType mean() const { return _mean; }
  RealType sigma() const { return _sigma; }

  void reset() { _valid = false; }

  template<class Engine>
  result_type operator()(Engine& eng)
  {
#ifndef BOOST_NO_STDC_NAMESPACE
    // allow for Koenig lookup
    using std::sqrt; using std::log;
#endif
    if(!_valid) {
      do
	{
	  _r1 = 2.0 * eng() - 1.0;
	  _r2 = 2.0 * eng() - 1.0;
	  _sq = _r1 *_r1 + _r2 *_r2;
	} while ((_sq > 1.0) || (_sq == 0.0));
      
      _sq = sqrt(-2 * log(1.0-_sq)/_sq);
      _r1 *= _sq;
      _r2 *= _sq;      
      _valid = true;
    } else
      _valid = false;
    
    return (_valid ? _r1 : _r2) * _sigma + _mean;
  }

#if !defined(BOOST_NO_OPERATORS_IN_NAMESPACE) && !defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)
  template<class CharT, class Traits>
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const normal_distribution& nd)
  {
    os << nd._mean << " " << nd._sigma << " "
       << nd._valid << " " << " " << nd._r1 << " " 
       << nd._r2;
    return os;
  }

  template<class CharT, class Traits>
  friend std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT,Traits>& is, normal_distribution& nd)
  {
    is >> std::ws >> nd._mean >> std::ws >> nd._sigma
       >> std::ws >> nd._valid >> std::ws >> nd._r1 
       >> std::ws >> nd._r2;
    return is;
  }
#endif
private:
  result_type _mean, _sigma;
  result_type _r1, _r2, _sq;
  bool _valid;
};

} // namespace boost

#endif // BOOST_RANDOM_NORMAL_DISTRIBUTION_HPP

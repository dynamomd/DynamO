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
#include <algorithm>

namespace magnet {
  namespace intersection {
    namespace detail {
      /*! \brief Implementation of unions of overlap functions. */
      template<class TA, class TB>
      class IntersectionOverlapFunction {
      public:
	IntersectionOverlapFunction(const TA& fA, const TB& fB):
	  _fA(fA), _fB(fB) {}
	
	inline double operator()(double dt = 0) const {
	  return std::min(_fA(dt), _fB(dt));
	}

	void timeShift(double dt) {
	  _fA.timeShift(dt);
	  _fB.timeShift(dt);
	}

	void flipSign() {
	  _fA.flipSign();
	  _fB.flipSign();
	}

	double nextEvent() const {
	  return std::min(_fA.nextEvent(), _fB.nextEvent());
	}

      private:
	TA _fA;
	TB _fB;
      };

      /*! \brief Helper function for creating unions of overlap
	functions. */
      template<class TA, class TB>
      IntersectionOverlapFunction<TA, TB> make_intersection(const TA& fA, const TB& fB)
      { return IntersectionOverlapFunction<TA, TB>(fA, fB); }
    }
  }
}

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
#include <magnet/math/quaternion.hpp>
#include <magnet/intersection/generic_algorithm.hpp>

namespace magnet {
  namespace intersection {
    namespace detail {
      /*! \brief The overlap function and its derivatives for two
	offcentre spheres rotating about an individual point.
      */
      class OffcentreSpheresOverlapFunction
      {
      public:
	OffcentreSpheresOverlapFunction(const math::Vector& rij, const math::Vector& vij, const math::Vector& omegai, const math::Vector& omegaj,
					const math::Vector& nu1, const math::Vector& nu2, const double diameter1, const double diameter2, 
					const double maxdist, const double t, const double invgamma, const double t_min, const double t_max):
	  w1(omegai), w2(omegaj), u1(nu1), u2(nu2), r12(rij), v12(vij), _diameter1(diameter1), _diameter2(diameter2), _invgamma(invgamma), _t(t), _t_min(t_min), _t_max(t_max)
	{
	  double Gmax = std::max(1 + t * invgamma, 1 + (t + t_max) * invgamma);
	  const double sigmaij = 0.5 * (_diameter1 + _diameter2);
	  const double sigmaij2 = sigmaij * sigmaij;
	  double magw1 = w1.nrm(), magw2 = w2.nrm();
	  double rijmax = Gmax * maxdist;
	  double magu1 = u1.nrm(), magu2 = u2.nrm();
	  double vijmax = v12.nrm() + Gmax * (magu1 * magw1 + magu2 * magw2) + std::abs(invgamma) * (magu1 + magu2);
	  double aijmax = Gmax * (magu1 * magw1 * magw1 + magu2 * magw2 * magw2) + 2 * std::abs(invgamma) * (magu1 * magw1 + magu2 * magw2);
	  double dotaijmax = Gmax * (magu1 * magw1 * magw1 * magw1 + magu2 * magw2 * magw2 * magw2) + 3 * std::abs(invgamma) * (magu1 * magw1 * magw1 + magu2 * magw2 * magw2);

	  _f1max = 2 * rijmax * vijmax + 2 * Gmax * std::abs(invgamma) * sigmaij2;
	  _f2max = 2 * vijmax * vijmax + 2 * rijmax * aijmax + 2 * invgamma * invgamma * sigmaij2;
	  _f3max = 6 * vijmax * aijmax + 2 * rijmax * dotaijmax;
	}

	template<size_t first_deriv=0, size_t nderivs = 1>
	std::array<double, nderivs> eval(const double dt = 0) const
	{
	  math::Vector u1new = Rodrigues(w1 * dt) * math::Vector(u1);
	  math::Vector u2new = Rodrigues(w2 * dt) * math::Vector(u2);

	  const double colldiam = 0.5 * (_diameter1 + _diameter2);
	  const double growthfactor = 1 + _invgamma * (_t + dt);
	  const math::Vector rij = r12 + dt * v12 + growthfactor * (u1new - u2new);
	  const math::Vector vij = v12 + growthfactor * ((w1 ^ u1new) - (w2 ^ u2new)) + _invgamma * (u1new - u2new);
	  const math::Vector aij = growthfactor * (-w1.nrm2() * u1new + w2.nrm2() * u2new) + 2 * _invgamma * ((w1 ^ u1new) - (w2 ^ u2new));
	  const math::Vector dotaij = growthfactor * (-w1.nrm2() * (w1 ^ u1new) + w2.nrm2() * (w2 ^ u2new)) + 3 * _invgamma * (-w1.nrm2() * u1new + w2.nrm2() * u2new);

	  std::array<double, nderivs> retval;
	  for (size_t i(0); i < nderivs; ++i)
	    switch (first_deriv + i) {
	    case 0: retval[i] = (rij | rij) - growthfactor * growthfactor * colldiam * colldiam; break;
	    case 1: retval[i] = 2 * (rij | vij) - 2 * _invgamma * growthfactor * colldiam * colldiam; break;
	    case 2: retval[i] = 2 * vij.nrm2() + 2 * (rij | aij) - 2 * _invgamma * _invgamma * colldiam * colldiam; break;
	    case 3: retval[i] = 6 * (vij | aij) + 2 * (rij | dotaij); break;
	    default:
	      M_throw() << "Invalid access";
	    }
	  return retval;
	}
    
	template<size_t deriv> 
	double max() const
	{
	  switch (deriv)
	    {
	    case 1: return _f1max;
	    case 2: return _f2max;
	    case 3: return _f3max;
	    default:
	      M_throw() << "Invalid access";
	    }
	}

	std::pair<bool, double> nextEvent() const {
	  return magnet::intersection::nextEvent(*this, _t_min, _t_max);
	}
  
      private:
	const math::Vector w1;
	const math::Vector w2;
	const math::Vector u1;
	const math::Vector u2;
	const math::Vector r12;
	const math::Vector v12;

	const double _diameter1, _diameter2, _invgamma;
	double _t, _f1max, _f2max, _f3max;
	const double _t_min, _t_max;
      };
    }
  }
}

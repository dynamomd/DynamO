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
					const double maxdist):
	  w1(omegai), w2(omegaj), u1(nu1), u2(nu2), r12(rij), v12(vij), _diameter1(diameter1), _diameter2(diameter2)
	{
	  double magw1 = w1.nrm(), magw2 = w2.nrm();
	  double rijmax = u1.nrm() + u2.nrm() + maxdist;
	  double vijmax = v12.nrm() + magw1 * u1.nrm() + magw2 * u2.nrm();
	  double aijmax = w1.nrm2() * u1.nrm() + w2.nrm2() * u2.nrm();
	  double dotaijmax = magw1 * w1.nrm2() * u1.nrm() + magw2 * w2.nrm2() * u2.nrm();
	  _f1max = 2 * rijmax * vijmax;
	  _f2max = 2 * vijmax * vijmax + 2 * rijmax * aijmax;
	  _f3max = 6 * vijmax * aijmax + 2 * rijmax * dotaijmax;
	}

	void stream(const double& dt)
	{
	  u1 = Rodrigues(w1 * dt) * math::Vector(u1);
	  u2 = Rodrigues(w2 * dt) * math::Vector(u2);
	  r12 += v12 * dt;
	}
  
	template<size_t deriv> 
	double eval() const
	{
	  double colldiam = 0.5 * (_diameter1 + _diameter2);
	  const math::Vector rij = r12 + u1 - u2;
	  const math::Vector vij = v12 + (w1 ^ u1) - (w2 ^ u2);
	  const math::Vector aij = -w1.nrm2() * u1 + w2.nrm2() * u2;
	  const math::Vector dotaij = -w1.nrm2() * (w1 ^ u1) + w2.nrm2() * (w2 ^ u2);

	  switch (deriv)
	    {
	    case 0: return (rij | rij) - colldiam * colldiam;
	    case 1: return 2 * (rij | vij);
	    case 2: return 2 * vij.nrm2() + 2 * (rij | aij);
	    case 3: return 6 * (vij | aij) + 2 * (rij | dotaij);
	    default:
	      M_throw() << "Invalid access";
	    }
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

	const math::Vector& getu1() const { return u1; }
	const math::Vector& getu2() const { return u2; }
	const math::Vector& getw1() const { return w1; }
	const math::Vector& getw2() const { return w2; }
	const math::Vector& getr12() const { return r12; }
	const math::Vector& getv12() const { return v12; }
  
	bool test_root() const { return true; }

      private:
	const math::Vector& w1;
	const math::Vector& w2;
	math::Vector u1;
	math::Vector u2;
	math::Vector r12;
	math::Vector v12;

	const double _diameter1, _diameter2;
	double _f1max, _f2max, _f3max;
      };

      /*! \brief The overlap function and its derivatives for two
	offcentre spheres rotating about an individual point.
      */
      class OffcentreGrowingSpheresOverlapFunction
      {
      public:
	OffcentreGrowingSpheresOverlapFunction(const math::Vector& rij, const math::Vector& vij, const math::Vector& omegai, const math::Vector& omegaj,
					       const math::Vector& nu1, const math::Vector& nu2, const double diameter1, const double diameter2, 
					       const double maxdist, const double t, const double invgamma, const double t_max):
	  w1(omegai), w2(omegaj), u1(nu1), u2(nu2), r12(rij), v12(vij), _diameter1(diameter1), _diameter2(diameter2), _invgamma(invgamma), _t_max(t_max), _t(t)
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

	void stream(const double& dt)
	{
	  u1 = Rodrigues(w1 * dt) * math::Vector(u1);
	  u2 = Rodrigues(w2 * dt) * math::Vector(u2);
	  _t += dt;
	  r12 += v12 * dt;
	}
  
	template<size_t deriv> 
	double eval() const
	{
	  const double colldiam = 0.5 * (_diameter1 + _diameter2);
	  const double growthfactor = 1 + _invgamma * _t;
	  const math::Vector rij = r12 + growthfactor * (u1 - u2);
	  const math::Vector vij = v12 + growthfactor * ((w1 ^ u1) - (w2 ^ u2)) + _invgamma * (u1 - u2);
	  const math::Vector aij = growthfactor * (-w1.nrm2() * u1 + w2.nrm2() * u2) + 2 * _invgamma * ((w1 ^ u1) - (w2 ^ u2));
	  const math::Vector dotaij = growthfactor * (-w1.nrm2() * (w1 ^ u1) + w2.nrm2() * (w2 ^ u2)) + 3 * _invgamma * (-w1.nrm2() * u1 + w2.nrm2() * u2);

	  switch (deriv)
	    {
	    case 0: return (rij | rij) - growthfactor * growthfactor * colldiam * colldiam;
	    case 1: return 2 * (rij | vij) - 2 * _invgamma * growthfactor * colldiam * colldiam;
	    case 2: return 2 * vij.nrm2() + 2 * (rij | aij) - 2 * _invgamma * _invgamma * colldiam * colldiam;
	    case 3: return 6 * (vij | aij) + 2 * (rij | dotaij);
	    default:
	      M_throw() << "Invalid access";
	    }
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

	const math::Vector& getu1() const { return u1; }
	const math::Vector& getu2() const { return u2; }
	const math::Vector& getw1() const { return w1; }
	const math::Vector& getw2() const { return w2; }
	const math::Vector& getr12() const { return r12; }
	const math::Vector& getv12() const { return v12; }
  
	bool test_root() const { return true; }

      private:
	const math::Vector& w1;
	const math::Vector& w2;
	math::Vector u1;
	math::Vector u2;
	math::Vector r12;
	math::Vector v12;

	const double _diameter1, _diameter2, _invgamma, _t_max;
	double _t, _f1max, _f2max, _f3max;
      };
    }

    /*! \brief Intersection test for offcentre spheres.
     */
    inline std::pair<bool, double> 
    offcentre_spheres(const math::Vector& rij, const math::Vector& vij, const math::Vector& angvi, const math::Vector& angvj,
		      const math::Vector& relativeposi, const math::Vector& relativeposj,
		      const double diameteri, const double diameterj, double maxdist, double t_max)
    {
#ifdef MAGNET_DEBUG
      if (std::isinf(t_max)) M_throw() << "Cannot perform root search in infinite intervals";
#endif

      detail::OffcentreSpheresOverlapFunction f(rij, vij, angvi, angvj, relativeposi, relativeposj, diameteri, diameterj, maxdist);
      return magnet::intersection::generic_algorithm(f, t_max, std::min(diameteri, diameterj) * 1e-10);
    }

    /*! \brief Intersection test for growing offcentre spheres.
     */
    inline std::pair<bool, double> 
    offcentre_growing_spheres(const math::Vector& rij, const math::Vector& vij, const math::Vector& angvi, const math::Vector& angvj,
			      const math::Vector& relativeposi, const math::Vector& relativeposj,
			      const double diameteri, const double diameterj, double maxdist, double t_max, double t, double invgamma)
    {
#ifdef MAGNET_DEBUG
      if (std::isinf(t_max)) M_throw() << "Cannot perform root search in infinite interval";
#endif

      detail::OffcentreGrowingSpheresOverlapFunction f(rij, vij, angvi, angvj, relativeposi, relativeposj, diameteri, diameterj, maxdist, t, invgamma, t_max);
      return magnet::intersection::generic_algorithm(f, t_max, std::min(diameteri, diameterj) * 1e-10);
    }
  }
}

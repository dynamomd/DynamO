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
#include <dynamo/particle.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/math/matrix.hpp>

namespace dynamo {
  template <class T, int derivative>
  class SFDerivative
  {
  public:
    SFDerivative(const T& sf): _shapeFunc(sf) {}

    void stream(const double& dt)
    { _shapeFunc.stream(dt); }

    template<size_t d>
    double eval() const
    { return _shapeFunc.template eval<d+derivative>(); }

    template<size_t d> 
    double max() const
    { return _shapeFunc.template max<d+derivative>(); }

    bool test_root() const { return true; }
    
  protected:
    T _shapeFunc;
  };

  /*! \brief The overlap function and its derivatives for offcentre
      spheres.
   */
  class SFOffcentre_Spheres
  {
  public:
    SFOffcentre_Spheres(const Vector& nr12, const Vector& nv12,
			const Vector& nw1, const Vector& nw2,
			const Quaternion& nq1, const Quaternion& nq2,
			const double diameter1, const double diameter2,
			const double offset1, const double offset2,
			const double maxdist):
      w1(nw1), w2(nw2), q1(nq1), q2(nq2),
      w12(nw1 - nw2), r12(nr12), v12(nv12),
      _diameter1(diameter1), _diameter2(diameter2),
      _offset1(offset1),_offset2(offset2)
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
      q1 = Quaternion::fromRotationAxis(w1 * dt) * q1;
      q1.normalise();
      q2 = Quaternion::fromRotationAxis(w2 * dt) * q2;
      q2.normalise();
      r12 += v12 * dt;
      u1 = q1 * Quaternion::initialDirector();
      u2 = q2 * Quaternion::initialDirector();
    }
  
    template<size_t deriv> 
    double eval() const
    {
      double colldiam = 0.5 * (_diameter1 + _diameter2);
      const Vector rij = r12 + u1 - u2;
      const Vector vij = v12 + (w1 ^ u1) - (w2 ^ u2);
      const Vector aij = -w1.nrm2() * u1 + w2.nrm2() * u2;
      const Vector dotaij = -w1.nrm2() * (w1 ^ u1) + w2.nrm2() * (w2 ^ u2);

      switch (deriv)
	{
	case 0:
	  return (rij | rij) - colldiam * colldiam;
	case 1:
	  return 2 * (rij | vij);
	case 2:
	  return 2 * vij.nrm2() + 2 * (rij | aij);
	case 3:
	  return 6 * (vij | aij) + 2 * (rij | dotaij);
	default:
	  M_throw() << "Invalid access";
	}
    }
    
    template<size_t deriv> 
    double max() const
    {
      switch (deriv)
	{
	case 1:
	  return _f1max;
	case 2:
	  return _f2max;
	case 3:
	  return _f3max;
	default:
	  M_throw() << "Invalid access";
	}
    }

    const Vector& getu1() const { return u1; }
    const Vector& getu2() const { return u2; }
    const Vector& getw1() const { return w1; }
    const Vector& getw2() const { return w2; }
    const Vector& getw12() const { return w12; }
    const Vector& getr12() const { return r12; }
    const Vector& getv12() const { return v12; }
  
    bool test_root() const { return true; }

  private:
    const Vector& w1;
    const Vector& w2;
    Quaternion q1;
    Quaternion q2;
    Vector u1;
    Vector u2;
    Vector w12;
    Vector r12;
    Vector v12;

    const double _diameter1, _diameter2;
    const double _offset1, _offset2;
    double _f1max, _f2max, _f3max;
  };
}

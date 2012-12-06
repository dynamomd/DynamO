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
#include <dynamo/dynamics/shapes/shape.hpp>
#include <magnet/math/matrix.hpp>

namespace dynamo {
  /*! \brief The overlap function and its derivatives for offcentre
      spheres.
      
      
   */
  class SFOffcentre_Spheres : public ShapeFunc {
  public:
    SFOffcentre_Spheres(const Vector& nr12, const Vector& nv12,
			const Vector& nw1, const Vector& nw2,
			const Vector& nu1, const Vector& nu2,
			const double diameter1, const double diameter2,
			const double maxdist):
      w1(nw1), w2(nw2), u1(nu1), u2(nu2),
      w12(nw1 - nw2), r12(nr12), v12(nv12),
      _diameter1(diameter1),
      _diameter2(diameter2)
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
      u1 = Rodrigues(w1 * dt) * Vector(u1);
      u2 = Rodrigues(w2 * dt) * Vector(u2);
      r12 += v12 * dt;
    }
  
    double F_zeroDeriv() const
    { 
      double colldiam = 0.5 * (_diameter1 + _diameter2);
      const Vector rij = r12 + u1 - u2;
      return (rij | rij) - colldiam * colldiam; 
    }

    double F_firstDeriv() const
    {    
      const Vector rij = r12 + u1 - u2;
      const Vector vij = v12 + (w1 ^ u1) - (w2 ^ u2);
      return 2 * (rij | vij);
    }

    double F_secondDeriv() const
    {
      const Vector rij = r12 + u1 - u2;
      const Vector vij = v12 + (w1 ^ u1) - (w2 ^ u2);
      const Vector aij = -w1.nrm2() * u1 + w2.nrm2() * u2;
      return 2 * vij.nrm2() + 2 * (rij | aij);
    }

    double F_thirdDeriv() const
    {
      const Vector rij = r12 + u1 - u2;
      const Vector vij = v12 + (w1 ^ u1) - (w2 ^ u2);
      const Vector aij = -w1.nrm2() * u1 + w2.nrm2() * u2;
      const Vector dotaij = -w1.nrm2() * (w1 ^ u1) + w2.nrm2() * (w2 ^ u2);
      return 6 * (vij | aij) + 2 * (rij | dotaij);
    }

    double F_firstDeriv_max() const { return _f1max; }
    double F_secondDeriv_max() const { return _f2max; }
    double F_thirdDeriv_max() const { return _f3max; }

    const Vector& getu1() const { return u1; }
    const Vector& getu2() const { return u2; }
    const Vector& getw1() const { return w1; }
    const Vector& getw2() const { return w2; }
    const Vector& getw12() const { return w12; }
    const Vector& getr12() const { return r12; }
    const Vector& getv12() const { return v12; }
  
    virtual bool test_root() const
    { return true; }

  private:
    const Vector& w1;
    const Vector& w2;
    Vector u1;
    Vector u2;
    Vector w12;
    Vector r12;
    Vector v12;

    const double _diameter1, _diameter2;
    double _f1max, _f2max, _f3max;
  };
}

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

#pragma once
#include "../../../simulation/particle.hpp"
#include "../../BC/BC.hpp"
#include "../../dynamics.hpp"
#include "../../../base/is_simdata.hpp"

class CLinesFunc {
public:
  CLinesFunc(const CParticle& np1, const CParticle& np2, 
	     const Vector& nw1, const Vector& nw2, 
	     const Vector& nu1, const Vector& nu2,
	     const DYNAMO::SimData* Sim):
    p1(np1), p2(np2), w1(nw1), w2(nw2), u1(nu1), u2(nu2),
    w12(nw1-nw2), r12(np1.getPosition() - np2.getPosition()),
    v12(np1.getVelocity() - np2.getVelocity())
  {
    Sim->Dynamics.BCs().setPBC(r12, v12);
  }

  void stream(const Iflt& dt)
  {
    u1 = Rodrigues(w1 * dt) * Vector(u1);
    u2 = Rodrigues(w2 * dt) * Vector(u2);
    r12 += v12 * dt;
  }
  
  std::pair<Iflt, Iflt> getCollisionPoints() const
  {
    Iflt rijdotui = (r12 | u1);
    Iflt rijdotuj = (r12 | u2);
    Iflt uidotuj = (u1 | u2);

    return std::make_pair(- (rijdotui - (rijdotuj * uidotuj)) / (1.0 - uidotuj*uidotuj),
			  (rijdotuj - (rijdotui * uidotuj)) / (1.0 - uidotuj*uidotuj));
  }

  Iflt F_zeroDeriv() const
  { return ((u1 ^ u2) | r12); }

  Iflt F_firstDeriv() const
  {    
    return ((u1 | r12) * (w12 | u2)) 
      + ((u2 | r12) * (w12 | u1)) 
      - ((w12 | r12) * (u1 | u2)) 
      + (((u1 ^ u2) | v12));
  }

  Iflt F_firstDeriv_max(const Iflt& length) const
  { return length * w12.nrm() + v12.nrm(); }

  Iflt F_secondDeriv(orientationStreamType A, orientationStreamType B) const
  {
    return 2.0 
      * (((u1 | v12) * (w12 | u2)) 
	 + ((u2 | v12) * (w12 | u1))
	 - ((u1 | u2) * (w12 | v12)))
      - ((w12 | r12) * (w12 | (u1 ^ u2))) 
      + ((u1 | r12) * (u2 | (w1 ^ w2))) 
      + ((u2 | r12) * (u1 | (w1 ^ w2)))
      + ((w12 | u1) * (r12 | (w2 ^ u2)))
      + ((w12 | u2) * (r12 | (w1 ^ u1))); 
  }

  Iflt F_secondDeriv_max(const Iflt& length) const
  {
    return w12.nrm() 
      * ((2 * v12.nrm()) + (length * (w1.nrm() + w2.nrm())));
  }

  const Vector& getu1() const { return u1; }
  const Vector& getu2() const { return u2; }
  const Vector& getw1() const { return w1; }
  const Vector& getw2() const { return w2; }
  const Vector& getw12() const { return w12; }
  const Vector& getr12() const { return r12; }
  const Vector& getv12() const { return v12; }
  
private:
  const CParticle& p1;
  const CParticle& p2;
  const Vector& w1;
  const Vector& w2;
  Vector u1;
  Vector u2;
  Vector w12;
  Vector r12;
  Vector v12;
};

/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "shape.hpp"

class COscillatingPlateFunc : public CShape {
public:
  COscillatingPlateFunc(const Vector& nvp, const Vector& nnhat,
		        const Vector& nrp, 
			const Iflt& nt, const Iflt& nDelta,
			const Iflt& nOmega, const Iflt& nSigma):
    vp(nvp), nhat(nnhat), rp(nrp),
    t(nt), Delta(nDelta), Omega(nOmega), Sigma(nSigma)
  {
    
  }

  void stream(const Iflt& dt)
  {
    t += dt;
    rp += vp * dt;
  }
  
  Iflt velnHatWall() const
  {
    Iflt retval = - Delta * Omega * std::sin(Omega * t);
    return  retval;
  }

  Iflt maxWallVel() const
  {    
    return  Delta * Omega;
  }

  Vector wallPosition() const
  {
    return  nhat * (Delta * std::cos(Omega * t));
  }

  Vector wallVelocity() const
  {
    return  -nhat * (Delta * Omega * std::sin(Omega * t));
  }

  Iflt F_zeroDeriv() const
  { 
    return ((rp - wallPosition()) | nhat) - Sigma;
  }

  Iflt F_zeroDerivFlip() const
  { 
    return ((rp - wallPosition()) | nhat) + Sigma;
  }

  Iflt F_firstDeriv() const
  {    
    return (vp | nhat) - velnHatWall();
  }

  Iflt F_firstDeriv_max(const Iflt&) const
  {
    return std::fabs(vp | nhat) + Delta * Omega; 
  }

  Iflt F_secondDeriv() const
  {
    return Delta * Omega * Omega * std::cos(Omega * t);
  }

  Iflt F_secondDeriv_max(const Iflt&) const
  {
    return Delta * Omega * Omega;
  }

  virtual CShape* Clone() const { return new COscillatingPlateFunc(*this); };

  virtual bool test_root(const Iflt&) const
  {
    return true;
  }

  void flipSigma() { Sigma = -Sigma; }
  
private:
  const Vector& vp;
  const Vector& nhat;
  Vector rp;
  lIflt t;
  Iflt Delta;
  Iflt Omega;
  Iflt Sigma;  
};

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

#ifndef CLNOrientation_H
#define CLNOrientation_H

#include "NewtonL.hpp"
#include <vector>
#include "../../datatypes/vector.hpp"

class CLinesFunc;

class CLNOrientation: public CLNewton
{
public:  
  CLNOrientation(DYNAMO::SimData* Sim, const XMLNode& XML):
    CLNewton(Sim),
    lastAbsoluteClock(-1),
    lastCollParticle1(0),
    lastCollParticle2(0)
  {}

  CLNOrientation(DYNAMO::SimData* Sim):
    CLNewton(Sim),
    lastAbsoluteClock(-1),
    lastCollParticle1(0),
    lastCollParticle2(0)
  {}

  virtual void initialise();

  virtual CLiouvillean* Clone() const { return new CLNOrientation(*this); }

  virtual void loadParticleXMLData(const XMLNode&, std::istream&);

  virtual void outputParticleBin64Data(std::ostream&) const;

  virtual void outputParticleXMLData(xmlw::XmlStream&) const;


  virtual bool getLineLineCollision(CPDData& PD, const Iflt& length, 
				    const CParticle& p1, const CParticle& p2
				    ) const;
  
  virtual C2ParticleData runLineLineCollision(const CIntEvent& eevent, 
					      const Iflt& elasticity, const Iflt& length) const;
  
  virtual C1ParticleData runAndersenWallCollision(const CParticle& part, 
						  const Vector & vNorm,
						  const Iflt& sqrtT
						  ) const;
  
  virtual C1ParticleData randomGaussianEvent(const CParticle& part, 
					     const Iflt& sqrtT) const;

  struct rotData
  {
    Vector  orientation;
    Vector  angularVelocity;
  };

  const rotData& getRotData(const CParticle& part) const
  { return orientationData[part.getID()]; }
  
  void initLineOrientations(const Iflt&);

protected:

  virtual void outputXML(xmlw::XmlStream&) const;

  //! \brief Helper function for writing out data
  template<class T>  void binarywrite(std::ostream&, const T&) const;
  template<class T>  void binaryread(std::istream&, T&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;
  
  virtual bool quadraticSolution(Iflt& returnVal, const int returnType, 
				 Iflt A, Iflt B, Iflt C) const;

  /* \brief For line line collisions, determines intersections of the infinite lines
  **
  **   Firstly, search for root in main window
  **  - If a root is not found, return failure
  **
  ** If a root is found: bring in an artificial new high boundary just beneath new root
  **  - If this leaves a window, search window for root
  **    - If a root is found, return to top of this section storing only this new root
  **    - If no root is found, drop out of this inner loop
  **  - Check root validity
  **    - If root is valid, this is earliest possible root - roll with it
  **    - If root is invalid, set new concrete t_low just above this found root and go from the top
  */
  virtual Iflt frenkelRootSearch(const CLinesFunc&, Iflt length, Iflt t_low, 
				 Iflt t_high) const;

  virtual Iflt quadraticRootHunter(const CLinesFunc& fL, Iflt length, 
				   Iflt& t_low, Iflt& t_high) const;

  virtual size_t getParticleDOF() const;
  virtual Iflt getParticleKineticEnergy(const CParticle& part) const;
  virtual void rescaleSystemKineticEnergy(const Iflt&);
  
  mutable std::vector<rotData> orientationData;
  mutable lIflt lastAbsoluteClock;
  mutable unsigned int lastCollParticle1;
  mutable unsigned int lastCollParticle2;  
};
#endif

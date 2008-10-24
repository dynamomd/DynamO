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

#ifndef CLiouvillean_H
#define CLiouvillean_H

#include "../../base/is_base.hpp"
#include "../eventtypes.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "datastruct.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CParticle;
class C2ParticleData;
class C1ParticleData;
class CIntEvent;
class intPart;

template <class T>
class CVector;
/*! \brief Describes the dynamics for single and pairs of objects
  between collisions.
 */
class CLiouvillean: public DYNAMO::SimBase
{
public:  
  CLiouvillean(DYNAMO::SimData* tmp):
    SimBase(tmp,"FreeStream", IC_blue),
    partPecTime(0.0),
    streamCount(0),
    streamFreq(1)
  {};

  virtual ~CLiouvillean() 
  {}

  inline void initialise() 
  {
    streamFreq = 10 * Sim->lN;
  }

  virtual bool SphereSphereInRoot(CPDData&, const Iflt&) const = 0;

  virtual bool SphereSphereOutRoot(CPDData&, const Iflt&) const = 0;  

  virtual bool sphereOverlap(const CPDData&, const Iflt&) const = 0;

  /*! \brief Determines when the particle center will hit a bounding box.

    Used by the cellular scheduler for cell transistion.
  */    
  virtual intPart getSquareCellCollision(const CParticle&, 
					 const CVector<>&, 
					 const CVector<>&) const = 0;

  /*! \brief Determines when the particle center will hit a bounding box.

    Used by the cellular global for cell transistion.
  */    
  virtual Iflt getSquareCellCollision2(const CParticle&, 
				       const CVector<>&, 
				       const CVector<>&) const = 0;
  
  /*! \brief Determines when the particle center will hit a bounding box.

    Used by the cellular scheduler for cell transistion.
  */    
  virtual size_t getSquareCellCollision3(const CParticle&, 
					 const CVector<>&, 
					 const CVector<>&) const = 0;

  /*! \brief Determines when the particle center will hit a wall.

    Used by the cellular scheduler for cell transistion. Will
    automatically wrap the box if travelling away from the nearest
    image of the wall.
  */    
  virtual Iflt getWallCollision(const CParticle&, 
				const CVector<>&, 
				const CVector<>&
				  ) const = 0;

  /*! \brief Determines when a particle has traveled half a box length
      in any dimension.

      This is used by CGSentinal to make sure no collisions are missed.
  */
  virtual Iflt getHalfBoxTraversalTime(const CParticle&) const = 0;

  /*! \brief Collides a particle with a wall.

    \param e Elasticity of wall
    \param vNorm Normal of the wall (\f$ vNorm \cdot v_1\f$ must be negative)
  */    
  virtual C1ParticleData runWallCollision(const CParticle&, 
					  const CVector<>& vNorm,
					  const Iflt& e
					  ) const = 0;

  /*! \brief Collides a particle with an Andersen thermostat wall.
    
    This gives a \f$ p \propto v_{norm} \exp(v_{norm}^2) \f$ distribution
    and gaussian tangent vectors

    \param sqrtT Square root of the Temperature of wall
    \param vNorm Normal of the wall (\f$ vNorm \cdot v_1 \f$ must be negative)
  */    
  virtual C1ParticleData runAndersenWallCollision(const CParticle&, 
						  const CVector<>& vNorm,
						  const Iflt& sqrtT
						  ) const = 0;
  
  /*! \brief Performs a hard sphere collision between the two particles.
    
    Also works for bounce collisions inside wells (i.e. will collide
    receeding particles).  
    \param e Elasticity
  */  
  virtual C2ParticleData SmoothSpheresColl(const CIntEvent&,const Iflt& e, const Iflt& d2, const EEventType& eType = CORE) const = 0;

  /*! \brief Does the maths for a well/shoulder event

      \param deltaKE kinetic energy change of event if executed
  */  
  virtual C2ParticleData SphereWellEvent(const CIntEvent&, const Iflt& deltaKE, const Iflt& d2) const = 0;

  virtual C1ParticleData randomGaussianEvent(const CParticle&, const Iflt&) const = 0;

  virtual CLiouvillean* Clone() const = 0;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CLiouvillean&);

  static CLiouvillean* loadClass(const XMLNode& ,DYNAMO::SimData*);
    
  void updateAllParticles() const
  {
    //May as well take this opportunity to reset the streaming
    //Note: the Replexing coordinator RELIES on this behaviour!
    BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
      {
	streamParticle(const_cast<CParticle&>(part), 
		       part.getPecTime() + partPecTime);
	
	const_cast<CParticle&>(part).getPecTime() = 0;
      }

    const_cast<Iflt&>(partPecTime) = 0;
    const_cast<size_t&>(streamCount) = 0;
  }

  inline void updateParticle(const CParticle& part) const
  {
    streamParticle(const_cast<CParticle&>(part), 
		   part.getPecTime() + partPecTime);

    const_cast<CParticle&>(part).getPecTime() = -partPecTime;
  }


  inline void updateParticlePair(const CParticle& p1, const CParticle& p2) const
  {
    //This is slow but sure, other stuff like reverse streaming, and
    //partial streaming are faster but work only for some collision
    //detections, not compression

    streamParticle(const_cast<CParticle&>(p1), 
		   partPecTime + p1.getPecTime());
    
    streamParticle(const_cast<CParticle&>(p2),
		   partPecTime + p2.getPecTime());
    
    const_cast<CParticle&>(p2).getPecTime()
      = const_cast<CParticle&>(p1).getPecTime() = -partPecTime;
  }

  inline void stream(const Iflt& dt)
  {
    partPecTime += dt;

    //Keep the magnitude of the partPecTime boundedx
    if (++streamCount == streamFreq)
      {
	BOOST_FOREACH(CParticle& part, Sim->vParticleList)
	  part.getPecTime() += partPecTime;

	partPecTime = 0;
	streamCount = 0;
      }
  }

  //This requires access to the streamParticle function and the
  //advanceUpdateParticle function
  friend class CSMultListShear;
  
protected:
  friend class CStreamTask;

  //See CSMultListShear, this just over advances the particle to find
  //its future position in boundary changes
  inline void advanceUpdateParticle(const CParticle& part, const Iflt& dt) const
  {
    streamParticle(const_cast<CParticle&>(part), dt + partPecTime + part.getPecTime());

    const_cast<CParticle&>(part).getPecTime() = - dt - partPecTime;
  }
  
  Iflt partPecTime;
  size_t streamCount;
  size_t streamFreq;
  
  virtual void outputXML(xmlw::XmlStream& ) const = 0;
  /*! \brief Moves the particle along in time */
  virtual void streamParticle(CParticle&, const Iflt&) const = 0;
};
#endif

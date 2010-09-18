/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CGNeighbourList_HPP
#define CGNeighbourList_HPP

#include "global.hpp"
#include <boost/function.hpp>
#include <vector>
#include <magnet/function/delegate.hpp>

class CGNeighbourList: public Global
{
public:
  typedef magnet::function::Delegate2
  <const Particle&, const size_t&, void> nbHoodFunc;

  typedef magnet::function::Delegate0<void> initFunc;  


protected:
  typedef std::pair<size_t, nbHoodFunc> nbHoodSlot;

  typedef std::pair<size_t, initFunc> initSlot;

  struct nbHoodSlotEraser
  {
    nbHoodSlotEraser(const size_t& id_in): 
      id(id_in)
    {}

    bool operator()(const nbHoodSlot& nbs) const
    { return  nbs.first == id; } 

    size_t id;
  };

  struct initSlotEraser
  {
    initSlotEraser(const size_t& id_in): 
      id(id_in)
    {}

    bool operator()(const initSlot& nbs) const
    { return  nbs.first == id; } 

    size_t id;
  };

public:
  CGNeighbourList(DYNAMO::SimData* a, 
		  const char *b): 
    Global(a,b)
  {}

  CGNeighbourList(CRange* a, DYNAMO::SimData* b, 
		  const char * c): 
    Global(a,b,c)
  {}

  CGNeighbourList(const CGNeighbourList&);

  virtual void getParticleNeighbourhood(const Particle&, 
					const nbHoodFunc&) const = 0;

  virtual void getParticleLocalNeighbourhood(const Particle&, 
					     const nbHoodFunc&) const = 0;

  template<class T> size_t
  ConnectSigCellChangeNotify
  (void (T::*func)(const Particle&, const size_t&)const , const T* tp) const 
  {    
    sigCellChangeNotify.push_back
      (nbHoodSlot(++sigCellChangeNotifyCount, 
		  getNBDelegate(func,tp)));
    
    return sigCellChangeNotifyCount; 
  }

  inline void
  DisconnectSigCellChangeNotify(const size_t& id) const 
  {    
    sigCellChangeNotify.erase
      (std::remove_if(sigCellChangeNotify.begin(),
		      sigCellChangeNotify.end(),
		      nbHoodSlotEraser(id)), 
      sigCellChangeNotify.end());
  }

  template<class T> size_t
  ConnectSigNewLocalNotify
  (void (T::*func)(const Particle&, const size_t&) const, const T* tp) const 
  {    
    sigNewLocalNotify.push_back
      (nbHoodSlot(++sigNewLocalNotifyCount, 
		  getNBDelegate(func,tp)));
    
    return sigNewLocalNotifyCount; 
  }

  inline void
  DisconnectSigNewLocalNotify(const size_t& id) const 
  {    
    sigNewLocalNotify.erase
      (std::remove_if(sigNewLocalNotify.begin(),
		      sigNewLocalNotify.end(),
		      nbHoodSlotEraser(id)), 
      sigNewLocalNotify.end());
  }


  template<class T> size_t
  ConnectSigNewNeighbourNotify
  (void (T::*func)(const Particle&, const size_t&) const, const T* tp) const 
  {    
    sigNewNeighbourNotify.push_back
      (nbHoodSlot(++sigNewNeighbourNotifyCount, 
		  getNBDelegate(func, tp)));
    
    return sigNewNeighbourNotifyCount; 
  }

  inline void
  DisconnectSigNewNeighbourNotify(const size_t& id) const 
  {    
    sigNewNeighbourNotify.erase
      (std::remove_if(sigNewNeighbourNotify.begin(),
		      sigNewNeighbourNotify.end(),
		      nbHoodSlotEraser(id)), 
      sigNewNeighbourNotify.end());
  }
    
  template<class T> size_t
  ConnectSigReInitNotify(void (T::*func)(), T* tp) const 
  {    
    sigReInitNotify.push_back
      (initSlot(++sigReInitNotifyCount, 
		getInitDelegate(func,tp)));
    
    return sigReInitNotifyCount; 
  }

  inline void
  DisconnectSigReInitNotify(const size_t& id) const 
  {    
    sigReInitNotify.erase
      (std::remove_if(sigReInitNotify.begin(),
		      sigReInitNotify.end(),
		      initSlotEraser(id)), 
      sigReInitNotify.end());
  }

  virtual Iflt getMaxSupportedInteractionLength() const = 0;

  virtual Iflt getMaxInteractionLength() const = 0;

  virtual void reinitialise(const Iflt&) = 0;

  void markAsUsedInScheduler() { isUsedInScheduler = true; }

protected:
  virtual void outputXML(xml::XmlStream&) const = 0;

  
  //Signals
  mutable size_t sigCellChangeNotifyCount;
  mutable std::vector<nbHoodSlot>
  sigCellChangeNotify;

  mutable size_t sigNewLocalNotifyCount;
  mutable std::vector<nbHoodSlot>
  sigNewLocalNotify;

  mutable size_t sigNewNeighbourNotifyCount;
  mutable std::vector<nbHoodSlot>
  sigNewNeighbourNotify;

  mutable size_t sigReInitNotifyCount;
  mutable std::vector<initSlot> 
  sigReInitNotify;

  bool isUsedInScheduler;
};

#endif

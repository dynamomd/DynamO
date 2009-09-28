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

#ifndef CSComplex_H
#define CSComplex_H

#include "scheduler.hpp"

class CSComplex: public CScheduler
{
public:
  CSComplex(const XMLNode&, DYNAMO::SimData* const);

  CSComplex(DYNAMO::SimData* const, CSSorter*);

  /*! \brief Must be overloaded to maintain connection status */
  CSComplex(const CSComplex&);

  virtual void rebuildList() { initialise(); }

  virtual void initialise();

  virtual void addEvents(const CParticle&);
  
  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addInteractionEvent(const CParticle&, const size_t&) const;
  void addInteractionEventInit(const CParticle&, const size_t&) const;

  void addEventsInit(const CParticle&);

  void addLocalEvent(const CParticle&, const size_t&) const;
  
  size_t NBListID;

  void virtualCellNewNeighbour(const CParticle&, const CParticle&);

  size_t cellChange;
  size_t cellChangeLocal;
  size_t reinit;         
};

#endif

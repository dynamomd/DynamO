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

#ifndef CSNeighbourList_H
#define CSNeighbourList_H

#include "scheduler.hpp"

class CSNeighbourList: public CScheduler
{
public:
  CSNeighbourList(const XMLNode&, DYNAMO::SimData* const);

  CSNeighbourList(DYNAMO::SimData* const, CSSorter*);

  /*! \brief Must be overloaded to maintain connection status */
  CSNeighbourList(const CSNeighbourList&);

  virtual void rebuildList() { initialise(); }

  virtual void initialise();

  virtual void addEvents(const CParticle&);
  
  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addEventsInit(const CParticle&);
  
  size_t NBListID;

  void virtualCellNewNeighbour(const CParticle&, const CParticle&);

  size_t cellChange;
  size_t cellChangeLocal;
  size_t reinit;         
};

#endif

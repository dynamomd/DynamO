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

#ifndef CSDumb_H
#define CSDumb_H

#include "scheduler.hpp"

class CSDumb: public CScheduler
{
public:
  CSDumb(const XMLNode&, const DYNAMO::SimData*);

  CSDumb(const DYNAMO::SimData*, CSSorter*);

  virtual void rebuildList() { initialise(); }

  virtual void initialise();

  virtual void update(const CParticle&);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addNewEvents(const CParticle&) const;
};

#endif

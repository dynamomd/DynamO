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

#ifndef CSMultListShear_H
#define CSMultListShear_H
#include "multlist.hpp"
#include "../datatypes/pluginpointer.hpp"
#include "../dynamics/ranges/1RRange.hpp"

class CSMultListShear: public CSMultList
{
public:
  CSMultListShear(const XMLNode&, const DYNAMO::SimData*);
  CSMultListShear(const DYNAMO::SimData*);

  void initialise();
  void reinitialise(Iflt);

  virtual void operator<<(const XMLNode&);

protected: 
  virtual void cellEvent(const CParticle&) const;
  virtual void outputXML(xmlw::XmlStream&) const;
};

#endif

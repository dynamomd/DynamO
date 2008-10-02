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

#ifndef CSMultListSpecial_H
#define CSMultListSpecial_H
#include "multlist.hpp"
#include "../datatypes/pluginpointer.hpp"
#include "../dynamics/ranges/1RRange.hpp"

class CSMultListSpecial: public CSMultList
{
public:
  CSMultListSpecial(const XMLNode&, const DYNAMO::SimData*);
  CSMultListSpecial(const DYNAMO::SimData*);

  void initialise();
  void reinitialise(Iflt);

  virtual void update(const CParticle&);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addSpecialEvents(const CParticle&);
  void addSpecialEvents_init(const CParticle&);

  CRRange specialParticles;
};

#endif

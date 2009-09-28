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

#ifndef CSCEntry_H
#define CSCEntry_H

#include "../../base/is_base.hpp"
#include "../../datatypes/pluginpointer.hpp"
#include "../../dynamics/ranges/1range.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}

class CParticle;

class CSCEntry: public DYNAMO::SimBase
{
public:
  CSCEntry(DYNAMO::SimData* const, const char *);
  
  virtual ~CSCEntry() = 0;

  virtual void initialise() = 0;
  
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CSCEntry&);

  static CSCEntry* getClass(const XMLNode&, DYNAMO::SimData* const);

  virtual void operator<<(const XMLNode&) = 0;
  
protected:
  
  smrtPlugPtr<CRange> range;
};

#endif

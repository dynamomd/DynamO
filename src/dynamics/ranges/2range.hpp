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

#ifndef C2Range_H
#define C2Range_H

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class Particle;

namespace DYNAMO
{
  class SimData;
}

class C2Range
{
public:
  virtual ~C2Range() {};
 
  virtual bool isInRange(const Particle&, const Particle&) const =0;  
  virtual void operator<<(const XMLNode& XML) = 0;
  
  virtual C2Range* Clone() const = 0;

  static C2Range* loadClass(const XMLNode&, const DYNAMO::SimData*);

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const C2Range&);

protected:
  virtual void outputXML(xmlw::XmlStream& XML) const = 0;
};

#endif

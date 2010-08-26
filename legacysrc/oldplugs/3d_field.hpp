/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef OP3DField_H
#define OP3DField_H

#include "outputplugin.hpp"
#include "../datatypes/field_array.hpp"

class OP3DField: public OutputPlugin
{
   public:
  OP3DField(DYNAMO::SimData*);
  ~OP3DField();

  void collisionUpdate(const IntEvent &, const CIntEventData &);
  
  void output(xml::XmlStream &);

  virtual OutputPlugin *Clone() const { return new OP3DField(*this); };
  
 protected:

  CFieldArray<Iflt> Density, Vsquared;
  CFieldArray<long> SampleCounter;
  CFieldArray<Vector  > Velocity;

  unsigned long imageCounter;
};

#endif

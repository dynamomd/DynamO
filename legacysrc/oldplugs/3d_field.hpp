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

#ifndef COP3DField_H
#define COP3DField_H

#include "outputplugin.hpp"
#include "../datatypes/field_array.hpp"

class COP3DField: public COutputPlugin
{
   public:
  COP3DField(DYNAMO::SimData*);
  ~COP3DField();

  void collisionUpdate(const CIntEvent &, const CIntEventData &);
  
  void output(xmlw::XmlStream &);

  virtual COutputPlugin *Clone() const { return new COP3DField(*this); };
  
 protected:

  CFieldArray<Iflt> Density, Vsquared;
  CFieldArray<long> SampleCounter;
  CFieldArray<Vector  > Velocity;

  unsigned long imageCounter;
};

#endif

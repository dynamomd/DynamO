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

#ifndef COPPRESSURE_H
#define COPPRESSURE_H

#include "outputplugin.hpp"
#include "../base/constants.hpp"

class COPPressure: public COutputPlugin
{
 public:
  COPPressure(DYNAMO::SimData*);
  ~COPPressure();

  void collisionUpdate(const CIntEvent &, const CIntEventData &b);

  void output(xmlw::XmlStream &); 

  void periodicOutput();

  virtual COutputPlugin *Clone() const 
    { return new COPPressure(*this); }
  
 protected:

  Iflt kpressure[NDIM][NDIM];
  Iflt cpressure[NDIM][NDIM];
  Iflt stream[NDIM][NDIM];
  Iflt t;
};

#endif

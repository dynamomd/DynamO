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

#ifndef COPCollVelDist_H
#define COPCollVelDist_H

#include "outputplugin.hpp"
#include "../datatypes/histogram.hpp"
#include "../base/constants.hpp"
#include "../dynamics/eventtypes.hpp"
#include <map>

class COPKE;

class CollVelHist: public C1DHistogram
{
public:
  CollVelHist():C1DHistogram(initVal) {}
  
  static Iflt initVal;
};

class COPCollVelDist: public COutputPlugin
{
 public:
  COPCollVelDist(DYNAMO::SimData*);

  void collisionUpdate(const CIntEvent &, const CIntEventData &b);

  void output(xmlw::XmlStream &); 

  virtual COutputPlugin *Clone() const 
    { return new COPCollVelDist(*this); }

  void setKEptr(COPKE *ptrKE2)
    {ptrKE = ptrKE2; };
  
 protected:
  std::map<EEventType, CollVelHist> vf;
  
  const COPKE* ptrKE; 
};

#endif

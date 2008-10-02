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

#ifndef COPRADDIST_H
#define COPRADDIST_H

#include "outputplugin.hpp"
#include <vector>

#define widthScale 0.025

class COPRadialDist: public COutputPlugin
{
 public:
  COPRadialDist(DYNAMO::SimData*);
  ~COPRadialDist();

  void collisionUpdate(const CIntEvent &, const CIntEventData &);

  void output(xmlw::XmlStream &);   

  virtual COutputPlugin *Clone() const
    { return new COPRadialDist(*this); }

 protected:

  long sampleCount;  

  std::vector<long> bin;

  Iflt binWidth;
  long binCount;
};

#endif

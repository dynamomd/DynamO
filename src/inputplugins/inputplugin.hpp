/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef INPUTPLUGIN_H
#define INPUTPLUGIN_H

#include "../base/is_base.hpp"
#include "../datatypes/vector.hpp"

class CInputPlugin: public DYNAMO::SimBase
{
 public:
  CInputPlugin(DYNAMO::SimData*, const char *aName, const char *aColor=IC_cyan);

  virtual ~CInputPlugin() {};

  virtual void initialise() {};

  //Rescaling system
  void rescaleVels(double val = 1.0);

  void zeroMomentum();  

  void setCOMVelocity(const Vector);

  void zeroCentreOfMass();  
  
  void setPackFrac(double);

  void mirrorDirection(unsigned int);

  void zeroVelComp(size_t);
  
 protected:  
};

#endif

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

#ifndef COPPovray_H
#define COPPovray_H

#include "ticker.hpp"

class COPPovray: public COPTicker
{
 public:
  COPPovray(const DYNAMO::SimData*);

  virtual COutputPlugin *Clone() const
  { return new COPPovray(*this); }

  virtual void initialise() {}

  virtual void stream(Iflt) {}

  virtual void ticker();
  
 protected:
  int frameCount;

  void printImage();
};

#endif

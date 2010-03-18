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

#ifndef OPTinkerXYZ_H
#define OPTinkerXYZ_H

#include "ticker.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>


class OPTinkerXYZ: public OPTicker
{
 public:
  OPTinkerXYZ(const DYNAMO::SimData*, const XMLNode&);
  ~OPTinkerXYZ();

  virtual OutputPlugin *Clone() const
  { return new OPTinkerXYZ(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();
  
  void operator<<(const XMLNode&);

 protected:
  int frameCount;
  bool fileOutput;
  bool liveOutput;
  bool blockForVMD;

  bool P1track;

  void *clientsock;
  void *sock;

  std::vector<float> coords;

  void printFileImage();
  void printLiveImage();
};

#endif

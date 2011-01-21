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

#ifndef OPStructureImaging_H
#define OPStructureImaging_H

#include "ticker.hpp"

class OPStructureImaging: public OPTicker
{
 public:
  OPStructureImaging(const DYNAMO::SimData*, const XMLNode&);

  virtual OutputPlugin *Clone() const
  { return new OPStructureImaging(*this); }

  virtual void initialise();

  virtual void stream(double) {}

  virtual void ticker();

  virtual void changeSystem(OutputPlugin*);

  virtual void operator<<(const XMLNode&);

  virtual void output(xml::XmlStream&);
  
 protected:

  void printImage();
  size_t id;
  size_t imageCount;
  std::vector<std::vector<Vector  > > imagelist;
  std::string structureName;
};

#endif

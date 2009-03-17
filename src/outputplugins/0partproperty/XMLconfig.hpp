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

#ifndef COPCONFIG_H
#define COPCONFIG_H

#include "../outputplugin.hpp"

class COPConfig: public COutputPlugin
{
 public:
  COPConfig(const DYNAMO::SimData*);
  ~COPConfig();
  
  virtual void initialise() {}

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&) {}

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&) {}

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&) {}

  virtual void eventUpdate(const CSystem&, const CNParticleData&, 
			   const Iflt&) {}

  virtual void output(xmlw::XmlStream &); 
  
  virtual COutputPlugin *Clone() const { return new COPConfig(*this); };

  void fileOutput(const char *);

  void setRounding();

  void setUncompressed() { compressedOutput = false; }

private:
  bool rounding;
  bool compressedOutput;
};

#endif

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

#ifndef OUTPUTPLUGIN_H
#define OUTPUTPLUGIN_H

#include "../base/is_base.hpp"
#include "../base/is_simdata.hpp"
#include <boost/foreach.hpp>
#include "../extcode/xmlParser.h"

class CIntEvent;
class CGlobEvent;
class XMLNode;
class C2ParticleData;
class C1ParticleData;
class CNParticleData;
class CSystem;

namespace xmlw
{
  class XmlStream;
}

class COutputPlugin: public DYNAMO::SimBase_const
{
public:
  COutputPlugin(const DYNAMO::SimData*, const char*, unsigned char order=100, const char *aColor=IC_blue);
  
  inline virtual ~COutputPlugin() {}
  
  virtual void initialise() = 0;
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&) = 0;
  
  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&) = 0;

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&) = 0;

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&) = 0;
  
  virtual COutputPlugin *Clone() const = 0;
  
  virtual void output(xmlw::XmlStream&);
  
  virtual void periodicOutput();
  
  static COutputPlugin* getPlugin(const XMLNode&, const DYNAMO::SimData*);
  static COutputPlugin* getPlugin(const std::string, const DYNAMO::SimData*);
  
  inline bool operator<(const COutputPlugin& OP) const
  { return updateOrder < OP.updateOrder; }
  
  inline bool operator>(const COutputPlugin& OP) const
  { return updateOrder > OP.updateOrder; }
  
  virtual void changeSystem(COutputPlugin*) 
  { D_throw() << "This plugin hasn't been prepared for changes of system\n Plugin " <<  name; }
  
  virtual void temperatureRescale(const Iflt&)
  {}
  
protected:
  DYNAMO::Colorise_Text_Stream_Operator I_Pcout() const;
  
  // This sets the order in which these things are updated
  // 0 is first
  // 100 is default
  // 250 is last
  //
  // Lets other plugins take data from plugins before/after they are updated
  unsigned char updateOrder;

private:

  template<class T> static COutputPlugin* 
  testGeneratePlugin(const DYNAMO::SimData*, const XMLNode&);
};

#endif

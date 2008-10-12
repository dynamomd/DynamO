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

#ifndef COPVACF_H
#define COPVACF_H

#include "../outputplugin.hpp"
#include <boost/circular_buffer.hpp>
#include <vector>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../dynamics/include.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../base/is_simdata.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"

class COPVACF: public COutputPlugin
{
 public:
  COPVACF(const DYNAMO::SimData*, const XMLNode&);

  virtual void operator<<(const XMLNode&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);
    
  virtual void output(xmlw::XmlStream&);

  virtual void initialise();

  virtual COutputPlugin *Clone() const { return new COPVACF(*this); }
   
protected:  
  
  virtual void newG(const C2ParticleData&);
  virtual void newG(const C1ParticleData&);
  virtual void newG(const CNParticleData&);
  
  virtual void accPass();

  Iflt getdt();

  std::vector<boost::circular_buffer<CVector<> > > G;
  std::vector<std::vector<CVector<> > > accG2;
  long count;
  Iflt dt, currentdt;
  
  size_t CorrelatorLength;
  size_t currCorrLen;
  bool notReady;
};

#endif

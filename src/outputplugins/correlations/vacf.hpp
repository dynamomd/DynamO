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

#ifndef OPVACF_H
#define OPVACF_H

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

class OPVACF: public OutputPlugin
{
 public:
  OPVACF(const DYNAMO::SimData*, const XMLNode&);

  virtual void operator<<(const XMLNode&);

  virtual void eventUpdate(const GlobalEvent&, const NEventData&);

  virtual void eventUpdate(const LocalEvent&, const NEventData&);

  virtual void eventUpdate(const System&, const NEventData&, const Iflt&);
  
  virtual void eventUpdate(const IntEvent&, const PairEventData&);
    
  virtual void output(xml::XmlStream&);

  virtual void initialise();

  virtual OutputPlugin *Clone() const { return new OPVACF(*this); }
   
protected:  
  
  virtual void newG(const PairEventData&);
  virtual void newG(const ParticleEventData&);
  virtual void newG(const NEventData&);
  
  virtual void accPass();

  Iflt getdt();

  std::vector<boost::circular_buffer<Vector  > > G;
  std::vector<std::vector<Vector  > > accG2;
  long count;
  Iflt dt, currentdt;
  
  size_t CorrelatorLength;
  size_t currCorrLen;
  bool notReady;
};

#endif

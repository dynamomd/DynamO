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

#ifndef OPSelfDiffusionOrientationalGK_H
#define OPSelfDiffusionOrientationalGK_H

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

class OPSelfDiffusionOrientationalGK: public OutputPlugin
{
 public:
  OPSelfDiffusionOrientationalGK(const DYNAMO::SimData*, const XMLNode&);

  virtual void operator<<(const XMLNode&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  virtual void eventUpdate(const IntEvent&, const C2ParticleData&);

  virtual void output(xmlw::XmlStream&);

  virtual void initialise();

  virtual OutputPlugin *Clone() const { return new OPSelfDiffusionOrientationalGK(*this); }

  typedef std::pair<Vector,Vector> VUpair;

protected:

  virtual void newG(const C2ParticleData&);
  virtual void newG(const C1ParticleData&);
  virtual void newG(const CNParticleData&);

  virtual void accPass();

  Iflt getdt();

  std::vector<boost::circular_buffer<VUpair> > G;
  std::vector<std::vector<Iflt> > accG2_parallel, accG2_perp;
  long count;
  Iflt dt, currentdt;

  size_t CorrelatorLength;
  size_t currCorrLen;
  bool notReady;
};

#endif

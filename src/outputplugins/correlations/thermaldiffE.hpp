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

#ifndef OPThermalDiffusionE_H
#define OPThermalDiffusionE_H

#include "../outputplugin.hpp"
#include "../../datatypes/vector.hpp"
#include <boost/circular_buffer.hpp>
#include "../0partproperty/misc.hpp"

/*! \brief The Correlator class for the Thermal Diffusion.*/
class OPThermalDiffusionE: public OutputPlugin
{
public:
  OPThermalDiffusionE(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

  virtual OutputPlugin* Clone() const 
  { return new OPThermalDiffusionE(*this); }

  virtual void eventUpdate(const GlobalEvent&, const NEventData&);
  virtual void eventUpdate(const LocalEvent&, const NEventData&);
  virtual void eventUpdate(const CSystem&, const NEventData&, const Iflt&);
  virtual void eventUpdate(const IntEvent&, const PairEventData&);

  void operator<<(const XMLNode&);
  
protected:
  boost::circular_buffer<Vector  > G;
  std::vector<Vector  > accG2;
  size_t count;
  Iflt dt, currentdt;
  Vector  constDelG, delG;
  size_t currlen;
  bool notReady;
  size_t CorrelatorLength;
  boost::circular_buffer<Vector  > Gsp1;
  Vector  constDelGsp1;
  Vector  delGsp1;
  size_t species1;
  Vector  sysMom;
  Iflt massFracSp1;

  std::string species1name;

  Iflt rescaleFactor();
  
  Vector  impulseDelG(const PairEventData&);
  Vector  impulseDelG(const NEventData&); 
  
  void updateConstDelG(const PairEventData&);
  void updateConstDelG(const ParticleEventData&);
  void updateConstDelG(const NEventData&);

  void stream(const Iflt);

  void newG();

  void accPass();
};

#endif


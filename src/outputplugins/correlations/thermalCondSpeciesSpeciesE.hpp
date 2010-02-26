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

#ifndef OPThermalConductivitySpeciesSpeciesE_H
#define OPThermalConductivitySpeciesSpeciesE_H

#include "../outputplugin.hpp"
#include "../../datatypes/vector.hpp"
#include <boost/circular_buffer.hpp>


/*! \brief The Correlator class for the Thermal Conductivity.*/
class OPThermalConductivitySpeciesSpeciesE: public OutputPlugin
{
public:
  OPThermalConductivitySpeciesSpeciesE(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

  virtual OutputPlugin* Clone() const 
  { return new OPThermalConductivitySpeciesSpeciesE(*this); }
  
  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&); 
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void operator<<(const XMLNode&);

protected:
  std::vector<boost::circular_buffer<Vector  > > G;
  std::vector<std::vector<Vector   > > accG2;
  size_t count;
  std::vector<Vector  > constDelG;
  std::vector<Vector  > delG;
  Iflt dt, currentdt;
  size_t currlen;
  bool notReady;

  size_t CorrelatorLength;

  void stream(const Iflt&);

  void newG();

  void accPass();

  Iflt rescaleFactor();
  
  void impulseDelG(const C2ParticleData&);

  void impulseDelG(const CNParticleData&);

  void updateConstDelG(const CNParticleData&);
  void updateConstDelG(const C2ParticleData&);
  void updateConstDelG(const C1ParticleData&);
};

#endif


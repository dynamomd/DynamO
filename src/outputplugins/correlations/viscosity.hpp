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

#ifndef COPVISCOSITY_H
#define COPVISCOSITY_H

#include "../outputplugin.hpp"
#include "../../datatypes/vector.hpp"
#include <boost/circular_buffer.hpp>

/*! \brief The Correlator class for the viscosity.*/
class COPViscosity: public COutputPlugin //COPCorrelator<CVector<CVector<> > >
{
public:
  COPViscosity(const DYNAMO::SimData*, const XMLNode& XML);

  virtual void initialise();

  virtual void output(xmlw::XmlStream &);
  
  virtual COutputPlugin* Clone() const { return new COPViscosity(*this); }
  
  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);
  
  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void stream(const Iflt& edt);

protected:
  inline CVector<CVector<> > impulseDelG(const C2ParticleData&);
  inline CVector<CVector<> > impulseDelG(const CNParticleData&);

  virtual void updateConstDelG(const C2ParticleData&);
  virtual void updateConstDelG(const C1ParticleData&);
  virtual void updateConstDelG(const CNParticleData&);
  
  void newG(CVector<CVector<> >);
  void accPass();

  CVector<CVector<> > avgTrace;
  boost::circular_buffer<CVector<CVector<> > > G;
  std::vector<CVector<CVector<> > > accG2;

  size_t count;
  Iflt dt, currentdt;
  CVector<CVector<> > constDelG, delG;

  size_t currlen;
  bool notReady;

  size_t CorrelatorLength;
};

#endif

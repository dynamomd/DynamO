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

#include "correlator.hpp"
#include "../../datatypes/vector.hpp"

/*! \brief The Correlator class for the viscosity.*/
class COPViscosity: public COPCorrelator<CVector<CVector<> > >
{
public:
  COPViscosity(const DYNAMO::SimData*, const XMLNode& XML);

  virtual void initialise();

  virtual void output(xmlw::XmlStream &);
  
  virtual COutputPlugin* Clone() const { return new COPViscosity(*this); }
  
protected:
  virtual Iflt rescaleFactor();
  
  inline CVector<CVector<> > impulseDelG(const C2ParticleData&);

  virtual void updateConstDelG(const C2ParticleData&);
  virtual void updateConstDelG(const C1ParticleData&);
  
  void newG(CVector<CVector<> >);

  CVector<CVector<> > avgTrace;
};

#endif

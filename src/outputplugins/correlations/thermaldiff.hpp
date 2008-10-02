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

#ifndef COPThermalDiffusion_H
#define COPThermalDiffusion_H

#include "correlator.hpp"
#include "../../datatypes/vector.hpp"

/*! \brief The Correlator class for the Thermal Conductivity.*/
class COPThermalDiffusion: public COPCorrelator<CVector<> >
{
public:
  COPThermalDiffusion(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

  virtual COutputPlugin* Clone() const { return new COPThermalDiffusion(*this); }
  
protected:
  virtual Iflt rescaleFactor();
  
  CVector<> impulseDelG(const C2ParticleData&);
  
  virtual void updateConstDelG(const C2ParticleData&);

  virtual void updateConstDelG(const C1ParticleData&);

  virtual void stream(const Iflt);

  virtual void newG();

  virtual void accPass();

  std::list<CVector<> > Gsp1;
  CVector<> constDelGsp1;
  CVector<> delGsp1;
  
  const CSpecies* species1;
  
  CVector<> sysMom;

  Iflt massFracSp1;

};

#endif


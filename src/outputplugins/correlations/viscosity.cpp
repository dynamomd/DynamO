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

#include "viscosity.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../datatypes/vector.xml.hpp"

COPViscosity::COPViscosity(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COutputPlugin(tmp,"Viscosity", 60),
  avgTrace(CVector<>(0.0)),
  count(0),
  dt(0),
  currentdt(0.0),
  constDelG(CVector<CVector<> >(CVector<>(0.0))), 
  delG(CVector<CVector<> >(CVector<>(0.0))),
  currlen(0),
  notReady(true),
  CorrelatorLength(100)
{}

void 
COPViscosity::initialise()
{
  Sim->getOutputPlugin<COPKEnergy>();
  Sim->getOutputPlugin<COPMisc>();

  G.resize(CorrelatorLength, CVector<CVector<> >(CVector<>(0.0)));

  accG2.resize(CorrelatorLength, CVector<CVector<> > (CVector<>(0.0)));

  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 50.0 / CorrelatorLength;
      else
	dt = 10.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->Dynamics.getkT()) * CorrelatorLength);
    }

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    constDelG += (part.getVelocity().dyad(part.getVelocity())) 
    * CVector<>(Sim->Dynamics.getSpecies(part).getMass());
}

void 
COPViscosity::eventUpdate(const CGlobEvent& iEvent, const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
COPViscosity::eventUpdate(const CSystem&, const CNParticleData& PDat, const Iflt& edt) 
{ 
  stream(edt);
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
COPViscosity::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& PDat)
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
COPViscosity::stream(const Iflt& edt)
{
  //Move the time forward
  //currentdt += edt;
  
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {
      delG += constDelG * CVector<>(dt - currentdt);
      newG (delG);
      currentdt += edt - dt;
      
      while (currentdt >= dt)
	{
	  delG = constDelG * CVector<>(dt);
	  currentdt -= dt;
	  newG(delG);
	}
      //Now calculate the start of the new delG
      delG = constDelG * CVector<>(currentdt);
    }
  else
    {
      currentdt += edt;
      delG += constDelG * CVector<>(edt);
    }
}

void 
COPViscosity::newG(CVector<CVector<> > Gval)
{
  avgTrace += Gval;

  G.push_front(Gval);

  if (notReady)
    {
      if (++currlen != CorrelatorLength)
	return;

      notReady = false;
    }

  accPass();
}


CVector<CVector<> >
COPViscosity::impulseDelG(const C2ParticleData& colldat)
{
  return colldat.particle1_.getDeltaP().dyad(colldat.rij);
}

CVector<CVector<> > 
COPViscosity::impulseDelG(const CNParticleData& ndat) 
{ 
  CVector<CVector<> > acc(CVector<>(0.0));
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    acc += impulseDelG(dat);
  
  return acc;
}

inline void 
COPViscosity::output(xmlw::XmlStream &XML)
{
  Iflt rescaleFactor = 1.0
    / (Sim->Dynamics.units().unitTime() 
       //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitViscosity() * 2.0 
       * Sim->getOutputPlugin<COPKEnergy>()->getAvgkT() 
       //Count has been taken out due to the extra averaging of the constant piece 
       * Sim->Dynamics.units().simVolume());
  
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt / Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size() / (Sim->getOutputPlugin<COPMisc>())->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor
      << xmlw::attr("SampleCount") << count
      << xmlw::attr("columns")
      << "t ";

  char name[3] = "xx";

  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      {
	name[0] = 'x' + i;
	name[1] = 'x' + j;
	XML << name << " ";       
      }
  
  CVector<CVector<> > traceAverage = avgTrace / CVector<>(((Iflt) G.size()) + ((Iflt) count));

  CVector<CVector<> > P = traceAverage / CVector<>(dt * Sim->Dynamics.units().simVolume());

  XML << xmlw::tag("Pressure")
      << P / CVector<>(Sim->Dynamics.units().unitPressure())
      << xmlw::endtag("Pressure");

  Iflt AvgPressure = 0.0;
  for (int iDim = 0; iDim < NDIM; iDim++)
    AvgPressure += P[iDim][iDim];
  
  XML << xmlw::tag("PressureVals")
      << xmlw::attr("AvgPressure")
      << AvgPressure / (NDIM * Sim->Dynamics.units().unitPressure())
    //Needs to use the KEnergy plugin
    /*      << xmlw::attr("AvgZ")
      << AvgPressure * Sim->Dynamics.units().simVolume() 
      / (NDIM * Sim->lN * getkT())*/
      << xmlw::endtag("PressureVals");
  
  XML << xmlw::chardata();
  
  for (unsigned int i = 0; i < accG2.size(); i++)
    {
      XML << (i+1) * dt / Sim->Dynamics.units().unitTime();
      for (int j = 0; j < NDIM; j++)
	for (int k = 0; k < NDIM; k++)
	  if (k==j)
	    XML << "\t" << ((accG2[i][j][k] / count) - pow(traceAverage[j][k]*(i+1),2)) * rescaleFactor ;
	  else
	    XML << "\t" << ((accG2[i][j][k] / count)) * rescaleFactor; 
      
      XML << "\n";
    }
  
  XML << xmlw::endtag("EinsteinCorrelator");
}

void 
COPViscosity::updateConstDelG(const C2ParticleData& PDat)
{
  //add the new value
  CVector<CVector<> > v1 = PDat.particle1_.getParticle().getVelocity()
    .dyad(PDat.particle1_.getParticle().getVelocity());
  CVector<CVector<> > v2 = PDat.particle2_.getParticle().getVelocity()
    .dyad(PDat.particle2_.getParticle().getVelocity());
  CVector<CVector<> > oldv1 = PDat.particle1_.getOldVel()
    .dyad(PDat.particle1_.getOldVel());
  CVector<CVector<> > oldv2 = PDat.particle2_.getOldVel()
    .dyad(PDat.particle2_.getOldVel());
  
  constDelG +=  ((v1 - oldv1) * CVector<>(PDat.particle1_.getSpecies().getMass()))
    + ((v2 - oldv2) * CVector<>(PDat.particle2_.getSpecies().getMass()));
}

void 
COPViscosity::updateConstDelG(const C1ParticleData& PDat)
{
  CVector<CVector<> > v1 = PDat.getParticle().getVelocity()
    .dyad(PDat.getParticle().getVelocity());
  CVector<CVector<> > oldv1 = PDat.getOldVel()
    .dyad(PDat.getOldVel());

  constDelG += ((v1 - oldv1) * CVector<>(PDat.getSpecies().getMass()));
}

void 
COPViscosity::updateConstDelG(const CNParticleData& ndat)
{
  BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
    updateConstDelG(dat);
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    updateConstDelG(dat);
}

void 
COPViscosity::accPass()
{
  ++count;
  CVector<CVector<> > sum(CVector<>(0.0));
  
  for (size_t i = 0; i < CorrelatorLength; ++i)
    {
      sum += G[i];
      accG2[i] += sum * sum;
    }
}

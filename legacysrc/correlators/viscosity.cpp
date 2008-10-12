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

COPViscosity::COPViscosity(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPCorrelator<CVector<CVector<> > >(tmp,"Viscosity", XML),
  avgTrace(0.0)
{
}

void 
COPViscosity::initialise()
{
  COPCorrelator<CVector<CVector<> > >::initialise();

  accG2.resize(CorrelatorLength, CVector<CVector<> > (0.0));
  dt = getdt();
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    constDelG += (part.getVelocity().dyad(part.getVelocity())) 
    * Sim->Dynamics.getSpecies(part).getMass();
}

Iflt 
COPViscosity::rescaleFactor() 
{
  return  1.0
    / (Sim->Dynamics.units().unitTime() //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitViscosity() * 2.0 
       * ptrEnergy->getAvgkT() //Count has been taken out due to the extra averaging of the constant piece 
       * Sim->Dynamics.units().simVolume());
}


void 
COPViscosity::newG(CVector<CVector<> > Gval)
{
  //This ensures the list stays at accumilator size
  if (G.size () == accG2.size ())
    {
      G.pop_back ();
      G.push_front (Gval);
      
      avgTrace += Gval;
      
      accPass ();
    }
  else
    {
      avgTrace += Gval;
    
      G.push_front (Gval);
      if (G.size () == accG2.size ())
	accPass();
    }
  
}


CVector<CVector<> >
COPViscosity::impulseDelG(const C2ParticleData& colldat)
{
  return colldat.particle1_.getDeltaP().dyad(colldat.rij);
}

inline void 
COPViscosity::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt / Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size() / ptrMisc->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor()
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
  
  CVector<CVector<> > traceAverage = avgTrace / (((Iflt) G.size()) + ((Iflt) count));

  CVector<CVector<> > P = traceAverage / (dt * Sim->Dynamics.units().simVolume());

  XML << xmlw::tag("Pressure")
      << P / Sim->Dynamics.units().unitPressure()
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
  
  Iflt Factor = rescaleFactor();

  XML << xmlw::chardata();
  
  for (unsigned int i = 0; i < accG2.size(); i++)
    {
      XML << (i+1) * dt / Sim->Dynamics.units().unitTime();
      for (int j = 0; j < NDIM; j++)
	for (int k = 0; k < NDIM; k++)
	  if (k==j)
	    XML << "\t" << ((accG2[i][j][k] / count) - pow(traceAverage[j][k]*(i+1),2)) * Factor ;
	  else
	    XML << "\t" << ((accG2[i][j][k] / count)) * Factor; 
      
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
  
  constDelG +=  ((v1 - oldv1) * PDat.particle1_.getSpecies().getMass())
    + ((v2 - oldv2) * PDat.particle2_.getSpecies().getMass());
}

void 
COPViscosity::updateConstDelG(const C1ParticleData& PDat)
{
  CVector<CVector<> > v1 = PDat.getParticle().getVelocity()
    .dyad(PDat.getParticle().getVelocity());
  CVector<CVector<> > oldv1 = PDat.getOldVel()
    .dyad(PDat.getOldVel());

  constDelG += ((v1 - oldv1) * PDat.getSpecies().getMass());
}

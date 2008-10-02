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

#include "pressure.hpp"
#include <boost/foreach.hpp>
#include "../dynamics/include.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../base/is_simdata.hpp"

COPPressure::COPPressure(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"Pressure"),
  t(0.0)
{  
  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      stream[i][j] = 0.0;

  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      {	
	cpressure[i][j] = 0.0;
	kpressure[i][j] = 0.0;
      }

  //Build the constant part of the stream
  BOOST_FOREACH( const CParticle & Part, Sim->vParticleList)
    for (int i = 0; i < NDIM; i++)
      for (int j = i; j < NDIM; j++)
	stream[i][j] += Sim->Dynamics.getSpecies(Part).getMass() 
	  * Part.getVelocity()[i] * Part.getVelocity()[j];
}

COPPressure::~COPPressure()
{}

void 
COPPressure::collisionUpdate(const CIntEvent &collision,
			     const CIntEventData &preColl)
{    
  Iflt dt = collision.getdt();
  t += dt;
  
  CVector<> veli = (collision.getParticle1()).getVelocity();
  CVector<> velj = (collision.getParticle2()).getVelocity();
  CVector<> rij = preColl.r12;
  
  CVector<> oldveli =  preColl.oldVelVec1;
  CVector<> oldvelj =  preColl.oldVelVec2;
  
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      {
	kpressure[i][j] += dt * (stream[i][j]);
	cpressure[i][j] += rij[j] * preColl.p1Species.getMass() 
	  * (veli[i] - oldveli[i]);	
      }
  
  //Update the streaming pressure
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      stream[i][j] += preColl.p1Species.getMass() * 
	(veli[i] * veli[j] - oldveli[i]*oldveli[j])
	+ preColl.p2Species.getMass() * 
	(velj[i] * velj[j] - oldvelj[i]*oldvelj[j]);
}
	
void
COPPressure::output(xmlw::XmlStream &XML)
{
  char name[4] = "Pxx";
  Iflt normFactor = Sim->Dynamics.units().unitLength()
    /(Sim->Dynamics.units().simVolume()*t);
  
  Iflt pressure[NDIM][NDIM];
  
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      pressure[i][j] = kpressure[i][j] + cpressure[i][j];
  
  XML << xmlw::tag("Pressure");
  XML << xmlw::tag("P") << xmlw::attr("val")
      << (pressure[0][0] + pressure [1][1] + pressure[2][2])*normFactor/3.0
      << xmlw::endtag("P");
  
  XML << xmlw::tag("FullTensor");
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      {
	name[1] = 'x' + i;
	name[2] = 'x' + j;
	XML << xmlw::tag(name) << xmlw::attr("val") 
	    << (pressure[i][j]) * normFactor 
	    << xmlw::endtag(name);
      }
  XML << xmlw::endtag("FullTensor");
  
  XML << xmlw::tag("KineticTensor");
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      {
	name[1] = 'x' + i;
	name[2] = 'x' + j;
	XML << xmlw::tag(name) << xmlw::attr("val") 
	    << (kpressure[i][j]) * normFactor 
	    << xmlw::endtag(name);
      }
  XML << xmlw::endtag("KineticTensor");


  XML << xmlw::endtag("Pressure")
      << xmlw::tag("NEMDvisc")
      << xmlw::tag("InPlaneVisc") << xmlw::attr("val")
      << (pressure[0][0] - pressure[1][1])*normFactor/2.0 
      << xmlw::endtag("InPlaneVisc")
      << xmlw::tag("OutPlaneVisc") << xmlw::attr("val")
      << (2.0*pressure[2][2] - (pressure[1][1] + pressure[0][0]))
          *normFactor/4.0 
      << xmlw::endtag("OutPlaneVisc") 
      << xmlw::endtag("NEMDvisc");
}

void
COPPressure::periodicOutput()
{}

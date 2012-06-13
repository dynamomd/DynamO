/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/systems/sysTicker.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace dynamo {

  Dynamics::Dynamics(dynamo::SimData* tmp): 
    SimBase(tmp, "Dynamics")
  {}

  void 
  Dynamics::addSystemTicker()
  {
    BOOST_FOREACH(shared_ptr<System>& ptr, Sim->systems)
      if (ptr->getName() == "SystemTicker")
	M_throw() << "System Ticker already exists";
  
    Sim->systems.push_back(shared_ptr<System>(new SysTicker(Sim, Sim->lastRunMFT, "SystemTicker")));
  }

  double
  Dynamics::calcInternalEnergy() const
  {
    double intECurrent = 0.0;

    BOOST_FOREACH(const shared_ptr<Interaction> & plugptr, 
		  Sim->interactions)
      intECurrent += plugptr->getInternalEnergy();

    return intECurrent;
  }

  double
  Dynamics::getSimVolume() const
  { 
    double vol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      vol *= Sim->primaryCellSize[iDim];
    return vol;
  }


  double
  Dynamics::getNumberDensity() const
  {
    return Sim->N / getSimVolume();
  }

  double 
  Dynamics::getPackingFraction() const
  {
    double volume = 0.0;
  
    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      BOOST_FOREACH(const size_t& ID, *(sp->getRange()))
      volume += sp->getIntPtr()->getExcludedVolume(ID);
  
    return  volume / getSimVolume();
  }

  void 
  Dynamics::setCOMVelocity(const Vector COMVelocity)
  {  
    Vector sumMV(0,0,0);
 
    long double sumMass(0);

    //Determine the discrepancy VECTOR
    BOOST_FOREACH(Particle & Part, Sim->particleList)
      {
	Vector  pos(Part.getPosition()), vel(Part.getVelocity());
	Sim->BCs->applyBC(pos,vel);
	double mass = Sim->species[Part].getMass(Part.getID());
	//Note we sum the negatives!
	sumMV -= vel * mass;
	sumMass += mass;
      }
  
    sumMV /= sumMass;
  
    sumMV += COMVelocity;

    BOOST_FOREACH(Particle & Part, Sim->particleList)
      Part.getVelocity() =  Part.getVelocity() + sumMV;
  }

  void 
  Dynamics::SystemOverlapTest()
  {
    Sim->liouvillean->updateAllParticles();

    std::vector<Particle>::const_iterator iPtr1, iPtr2;
  
    for (iPtr1 = Sim->particleList.begin(); iPtr1 != Sim->particleList.end(); ++iPtr1)
      for (iPtr2 = iPtr1 + 1; iPtr2 != Sim->particleList.end(); ++iPtr2)    
	Sim->getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      BOOST_FOREACH(const shared_ptr<Local>& lcl, Sim->locals)
      if (lcl->isInteraction(part))
	lcl->checkOverlaps(part);
  }
}

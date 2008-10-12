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

#ifndef COPCORRELATOR_H
#define COPCORRELATOR_H

#include "../outputplugin.hpp"
#include <list>
#include <vector>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../dynamics/include.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../base/is_simdata.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"

template<typename T = Iflt>
class COPCorrelator: public COutputPlugin
{
 public:
  COPCorrelator(const DYNAMO::SimData* tmp,const char *aName, const XMLNode& XML):
    COutputPlugin(tmp, aName, 60), //Note the sort order set later
    count(0),
    dt(0),
    currentdt(0.0),
    constDelG(T(0.0)), 
    delG(T(0.0)),
    ptrEnergy(NULL),
    ptrMisc(NULL),
    CorrelatorLength(100)
  {
    operator<<(XML);
  }

  virtual void operator<<(const XMLNode& XML)
  {
    try 
      {
	if (XML.isAttributeSet("Length"))
	  CorrelatorLength = boost::lexical_cast<unsigned int>(XML.getAttribute("Length"));
	
	if (XML.isAttributeSet("dt"))
	  dt = Sim->Dynamics.units().unitTime() * 
	    boost::lexical_cast<Iflt>(XML.getAttribute("dt"));

	if (XML.isAttributeSet("t"))
	  dt = Sim->Dynamics.units().unitTime() * 
	    boost::lexical_cast<Iflt>(XML.getAttribute("t"))/CorrelatorLength;
      }
    catch (boost::bad_lexical_cast &)
      {
	D_throw() << "Failed a lexical cast in COPCorrelator";
      }

  }

  virtual void stream(const Iflt edt)
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

  virtual void eventUpdate(const CGlobEvent& iEvent, const CNParticleData& PDat) 
  {
    stream(iEvent.getdt());
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  virtual void eventUpdate(const CSystem&, const CNParticleData& PDat, const Iflt& edt) 
  { 
    stream(edt);
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }
  
  virtual void eventUpdate(const CIntEvent& iEvent, const C2ParticleData& PDat)
  {
    stream(iEvent.getdt());
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }
    
  virtual void output(xmlw::XmlStream& XML)
  {
    XML << xmlw::tag("Correlator")
	<< xmlw::attr("name") << name
	<< xmlw::attr("size") << accG2.size()
	<< xmlw::attr("dt") << dt / Sim->Dynamics.units().unitTime()
	<< xmlw::attr("LengthInMFT") << dt * accG2.size() / ptrMisc->getMFT()
	<< xmlw::attr("simFactor") << rescaleFactor()
	<< xmlw::attr("SampleCount") << count;
    
    Iflt factor = rescaleFactor();

    for (unsigned int i = 0; i < accG2.size(); i++)
      XML << xmlw::tag("data") << xmlw::attr("t")
	  << (i+1) * dt / Sim->Dynamics.units().unitTime()
	  << accG2[i] * CVector<>(factor)
	  << xmlw::endtag("data");
    
    XML << xmlw::endtag("Correlator");
  }

  virtual void initialise()
  {
    ptrEnergy = Sim->getOutputPlugin<COPKEnergy>();
    ptrMisc = Sim->getOutputPlugin<COPMisc>();
  }
 
  std::list<T> getAvgAcc() const
  {
    std::list<T> tmp;
    
    BOOST_FOREACH(const T &val, accG2)
      tmp.push_back(val/((Iflt) count));
    
    return tmp;
  }
  
 protected:  
  virtual T impulseDelG(const C2ParticleData&) { return T(0.0); }

  virtual T impulseDelG(const C1ParticleData&) { return T(0.0); }
  
  virtual T impulseDelG(const CNParticleData& ndat) 
  { 
    T acc(0);
    
    BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
      acc += impulseDelG(dat);
    
    BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
      acc += impulseDelG(dat);

    return acc;
  }

  virtual void updateConstDelG(const C2ParticleData&) {}

  virtual void updateConstDelG(const C1ParticleData&) {}

  virtual void updateConstDelG(const CNParticleData& ndat)
    {
      BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
	updateConstDelG(dat);

      BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
	updateConstDelG(dat);
    }

  virtual Iflt rescaleFactor() = 0;
  
  virtual void newG(T Gval)
  {
    //This ensures the list stays at accumilator size
    if (G.size () == accG2.size ())
      {
	G.pop_back ();
	G.push_front (Gval);
	accPass ();
      }
    else
      {
	G.push_front (Gval);
	if (G.size () == accG2.size ())
	  accPass();
      }
  }
  
  virtual void accPass()
  {
    count++;
    T sum(0);
    
    int index = 0;
    
    BOOST_FOREACH(const T &val, G) 
      {
	sum += val;
	accG2[index] += sum * sum;
	index++;
      }
  }

  Iflt getdt()
  {
    //Get the simulation temperature
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  return Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  return 10.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->Dynamics.getkT()) * CorrelatorLength);
      }
    else 
      return dt;
  }
    
  std::list<T> G;
  std::vector<T> accG2;
  long count;
  Iflt dt, currentdt;
  T constDelG, delG;

  const COPKEnergy* ptrEnergy;
  const COPMisc* ptrMisc;
  unsigned int CorrelatorLength;
};

#endif

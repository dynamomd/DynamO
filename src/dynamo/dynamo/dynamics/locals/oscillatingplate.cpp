/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/Function.hpp>
#endif 

#include "oscillatingplate.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLOscillatingPlate::CLOscillatingPlate(dynamo::SimData* nSim,
				       Vector nrw0, Vector nnhat,
				       double nomega0, double nsigma, double ne,
				       double ndelta, double nmass, std::string nname, 
				       CRange* nRange, double timeshift, bool nstrongPlate):
  Local(nRange, nSim, "OscillatingPlate"),
  strongPlate(nstrongPlate),
  rw0(nrw0), nhat(nnhat), omega0(nomega0), sigma(nsigma), 
  e(ne), delta(ndelta), mass(nmass), timeshift(0), 
  lastID(std::numeric_limits<size_t>::max()), lastdSysTime(HUGE_VAL)
{
  localName = nname;
}

CLOscillatingPlate::CLOscillatingPlate(const magnet::xml::Node& XML, dynamo::SimData* tmp):
  Local(tmp, "OscillatingPlate"),
  lastID(std::numeric_limits<size_t>::max()), lastdSysTime(HUGE_VAL)
{
  operator<<(XML);
}

LocalEvent 
CLOscillatingPlate::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  //bool caution = ((part.getID() == lastID) && (lastdSysTime == Sim->dSysTime));

  double reducedt = Sim->dSysTime 
    - 2.0 * M_PI * int(Sim->dSysTime * omega0 / (2.0*M_PI)) / omega0;

  std::pair<bool, double> eventData = Sim->dynamics.getLiouvillean().getPointPlateCollision
    (part, rw0, nhat, delta, omega0, sigma, reducedt + timeshift, 
     false);

  EEventType type = (eventData.first) ? WALL : VIRTUAL ;

  if (eventData.second == HUGE_VAL)
    type = NONE;
  
  return LocalEvent(part, eventData.second, type, *this);
}



void
CLOscillatingPlate::runEvent(const Particle& part, const LocalEvent& iEvent) const
{
  ++Sim->eventCount;
  
  //Run the collision and catch the data
  NEventData EDat(Sim->dynamics.getLiouvillean().runOscilatingPlate
		      (part, rw0, nhat, delta, omega0, sigma, mass, 
		       e, timeshift, strongPlate));

  lastdSysTime = Sim->dSysTime;
  lastID = part.getID();

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event update the scheduler and plugins
  //if (strongPlate) 
  //  Sim->ptrScheduler->fullUpdate(part);
  //else
    Sim->ptrScheduler->rebuildList();

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

bool 
CLOscillatingPlate::isInCell(const Vector & Origin, const Vector& CellDim) const
{
  return true;
}

void 
CLOscillatingPlate::initialise(size_t nID)
{
  ID = nID;
}

void 
CLOscillatingPlate::operator<<(const magnet::xml::Node& XML)
{
  range.set_ptr(CRange::getClass(XML,Sim));
  
  try {
    e = XML.getAttribute("Elasticity").as<double>();
    nhat << XML.getNode("Norm");
    nhat /= nhat.nrm();

    rw0 << XML.getNode("Origin");
    rw0 *= Sim->dynamics.units().unitLength();

    if (XML.getAttribute("StrongPlate").valid())
      strongPlate = XML.getAttribute("StrongPlate").as<double>();

    omega0 = XML.getAttribute("Omega0").as<double>() / Sim->dynamics.units().unitTime();
    sigma = XML.getAttribute("Sigma").as<double>() * Sim->dynamics.units().unitLength();
    delta = XML.getAttribute("Delta").as<double>() * Sim->dynamics.units().unitLength();
    mass = XML.getAttribute("Mass").as<double>()  * Sim->dynamics.units().unitMass();
    timeshift = XML.getAttribute("TimeShift").as<double>() * Sim->dynamics.units().unitTime();

    localName = XML.getAttribute("Name");
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CLOscillatingPlate";
    }
}

void 
CLOscillatingPlate::outputXML(xml::XmlStream& XML) const
{
  double tmp = Sim->dSysTime + timeshift;

  tmp -= 2.0 * M_PI * int(tmp * omega0 / (2.0 * M_PI) ) / omega0;

  XML << xml::attr("Type") << "OscillatingPlate" 
      << xml::attr("Name") << localName
      << xml::attr("Elasticity") << e
      << xml::attr("Omega0") << omega0 * Sim->dynamics.units().unitTime()
      << xml::attr("Sigma") << sigma / Sim->dynamics.units().unitLength()
      << xml::attr("Delta") << delta / Sim->dynamics.units().unitLength()
      << xml::attr("Mass") << mass / Sim->dynamics.units().unitMass()
      << xml::attr("TimeShift") << tmp / Sim->dynamics.units().unitTime()
      << xml::attr("StrongPlate") << strongPlate
      << range
      << xml::tag("Norm")
      << nhat
      << xml::endtag("Norm")
      << xml::tag("Origin")
      << rw0 / Sim->dynamics.units().unitLength()
      << xml::endtag("Origin");

}

Vector
CLOscillatingPlate::getPosition() const
{
  return nhat * (delta * std::cos(omega0 * (Sim->dSysTime + timeshift))) + rw0;
}

Vector
CLOscillatingPlate::getVelocity() const
{
  return - nhat * (delta * omega0 * std::sin(omega0 * (Sim->dSysTime + timeshift)));
}

double 
CLOscillatingPlate::getPlateEnergy() const
{
  return 0.5 * mass 
    * (std::pow(omega0 * delta * std::cos(omega0 * (Sim->dSysTime + timeshift)), 2)
       +  std::pow(omega0 * delta * std::sin(omega0 * (Sim->dSysTime + timeshift)), 2));
}

#ifdef DYNAMO_visualizer

magnet::thread::RefPtr<RenderObj>& 
CLOscillatingPlate::getCoilRenderObj() const
{
  const double lengthRescale = 1 / Sim->primaryCellSize.maxElement();

  if (!_renderObj.isValid())
    {
      Vector axis3 = nhat / nhat.nrm();
      Vector axis2(0,0,1);
      
      for (size_t i(0); i < NDIM; ++i)
	{
	  Vector tryaxis = Vector(0,0,0);
	  tryaxis[i] = 1;
	  Vector tryaxis2 = axis3 ^ tryaxis;
	  
	  if (tryaxis2.nrm() != 0) { axis2 = tryaxis2 / tryaxis2.nrm(); break; }
	}

      Vector axis1 = axis2 ^ axis3;

      std::ostringstream os;
      os << axis3[0] << ", "
	 << axis3[1] << ", "
	 << axis3[2] << ", 0";

      axis1 *= Sim->primaryCellSize[1] * lengthRescale / axis1.nrm();
      axis2 *= Sim->primaryCellSize[2] * lengthRescale / axis2.nrm();

      _renderObj = new RFunction(10, 
				 rw0 - 0.5 * (axis1 + axis2), 
				 axis1, axis2, axis3,
				 0, 0, 1, 1, true, false,
				 "Oscillating wall",
				 "f = A;",
				 "normal = -(float4)(" + os.str() + ");"
				 );
    }
  
  return _renderObj;
}

void 
CLOscillatingPlate::updateRenderData(magnet::CL::CLGLState&) const
{
  const double lengthRescale = 1 / Sim->primaryCellSize.maxElement();

  if (_renderObj.isValid())
    static_cast<RFunction&>(*_renderObj)
      .setConstantA((delta * std::cos(omega0 * (Sim->dSysTime + timeshift)) 
		     - (sigma + 0.5 * Sim->dynamics.units().unitLength())) *  lengthRescale);
}
#endif

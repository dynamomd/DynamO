/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "oscillatingplate.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLOscillatingPlate::CLOscillatingPlate(DYNAMO::SimData* nSim,
				       Vector nrw0, Vector nnhat,
				       double nomega0, double nsigma, double ne,
				       double ndelta, double nmass, std::string nname, 
				       CRange* nRange, double timeshift, bool nstrongPlate):
  Local(nRange, nSim, "OscillatingPlate"),
  strongPlate(nstrongPlate),
  rw0(nrw0), nhat(nnhat), omega0(nomega0), sigma(nsigma), 
  e(ne), delta(ndelta), mass(nmass), timeshift(0), lastID(-1), lastdSysTime(HUGE_VAL)
{
  localName = nname;
}

CLOscillatingPlate::CLOscillatingPlate(const XMLNode& XML, DYNAMO::SimData* tmp):
  Local(tmp, "OscillatingPlate"),
  lastID(-1), lastdSysTime(HUGE_VAL)
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
CLOscillatingPlate::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<double>(XML.getAttribute("Elasticity"));

    XMLNode xBrowseNode = XML.getChildNode("Norm");
    nhat << xBrowseNode;
    nhat /= nhat.nrm();

    xBrowseNode = XML.getChildNode("Origin");
    rw0 << xBrowseNode;
    rw0 *= Sim->dynamics.units().unitLength();

    if (XML.isAttributeSet("StrongPlate"))
      strongPlate = boost::lexical_cast<double>(XML.getAttribute("StrongPlate"));

    omega0 =  boost::lexical_cast<double>(XML.getAttribute("Omega0"))
      / Sim->dynamics.units().unitTime();

    sigma = boost::lexical_cast<double>(XML.getAttribute("Sigma"))
      * Sim->dynamics.units().unitLength();

    delta = boost::lexical_cast<double>(XML.getAttribute("Delta"))
      * Sim->dynamics.units().unitLength();

    mass = boost::lexical_cast<double>(XML.getAttribute("Mass"))
      * Sim->dynamics.units().unitMass();

    timeshift = boost::lexical_cast<double>(XML.getAttribute("TimeShift"))
      * Sim->dynamics.units().unitTime();

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

void 
CLOscillatingPlate::write_povray_info(std::ostream& os) const
{
  Vector pos = getPosition();
  //The walls are \pm0.25 thick and are set 0.5 away from the COM surface
  //to give the appearance of proper walls, not COM walls
  Vector WallLoc1 = pos + nhat * (sigma + 0.75 * Sim->dynamics.units().unitLength());
  Vector WallLoc2 = pos - nhat * (sigma + 0.75 * Sim->dynamics.units().unitLength());

  Sim->dynamics.BCs().applyBC(WallLoc1);
  Sim->dynamics.BCs().applyBC(WallLoc2);

  os << "object { intersection { object { box { <-0.5, " << -0.25 * Sim->dynamics.units().unitLength() 
     << ", -0.5>, <0.5, " << +0.25 * Sim->dynamics.units().unitLength() 
     << ", 0.5> } Point_At_Trans(<"
     << nhat[0] << "," << nhat[1] << "," << nhat[2] << ">) translate <"
     <<  WallLoc1[0] << "," <<  WallLoc1[1] << "," <<  WallLoc1[2] 
     << "> }"

     << "\nbox { <" 
     << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
     << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
     << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
     << ">,"
     << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
     << "> }\n"
     << "} pigment { Col_Glass_Bluish } finish { F_Glass5 } }\n";

  os << "object { intersection { object { box { <-0.5, " << -0.25 * Sim->dynamics.units().unitLength()
     << ", -0.5>, <0.5, " << 0.25 * Sim->dynamics.units().unitLength() 
     << ", 0.5> } Point_At_Trans(<"
     << -nhat[0] << "," << -nhat[1] << "," << -nhat[2] << ">) translate <"
     <<  WallLoc2[0] << "," <<  WallLoc2[1] << "," <<  WallLoc2[2] 
     << "> }"

     << "\nbox { <" 
     << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
     << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
     << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
     << ">,"
     << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
     << "> }\n"
     << "} pigment { Col_Glass_Bluish } finish { F_Glass5 } }\n";
}

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/Function.hpp>

magnet::thread::RefPtr<RenderObj>& 
CLOscillatingPlate::getCoilRenderObj() const
{
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

      axis1 *= Sim->aspectRatio[1] / axis1.nrm();
      axis2 *= Sim->aspectRatio[2] / axis2.nrm();

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
  if (_renderObj.isValid())
    static_cast<RFunction&>(*_renderObj)
      .setConstantA(delta * std::cos(omega0 * (Sim->dSysTime + timeshift)) 
		    - (sigma + 0.5 * Sim->dynamics.units().unitLength()));
}
#endif

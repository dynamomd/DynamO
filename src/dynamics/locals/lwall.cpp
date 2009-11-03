/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "lwall.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLWall::CLWall(DYNAMO::SimData* nSim, Iflt ne, Vector  nnorm, 
	       Vector  norigin, std::string nname, 
	       CRange* nRange, bool nrender):
  CLocal(nRange, nSim, "LocalWall"),
  vNorm(nnorm),
  vPosition(norigin),
  e(ne),
  render(nrender)
{
  localName = nname;
}

CLWall::CLWall(const XMLNode& XML, DYNAMO::SimData* tmp):
  CLocal(tmp, "LocalWall")
{
  operator<<(XML);
}

CLocalEvent 
CLWall::getEvent(const CParticle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  return CLocalEvent(part, Sim->dynamics.getLiouvillean().getWallCollision
		     (part, vPosition, vNorm), WALL, *this);
}

void
CLWall::runEvent(const CParticle& part, const CLocalEvent& iEvent) const
{
  ++Sim->lNColl;

  //Run the collision and catch the data
  CNParticleData EDat(Sim->dynamics.getLiouvillean().runWallCollision
		      (part, vNorm, e));

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

bool 
CLWall::isInCell(const Vector & Origin, const Vector& CellDim) const
{
  return DYNAMO::OverlapFunctions::CubePlane
    (Origin, CellDim, vPosition, vNorm);
}

void 
CLWall::initialise(size_t nID)
{
  ID = nID;
}

void 
CLWall::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
    render = boost::lexical_cast<bool>(XML.getAttribute("Render"));
    XMLNode xBrowseNode = XML.getChildNode("Norm");
    localName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm /= vNorm.nrm();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->dynamics.units().unitLength();
  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CLWall";
    }
}

void 
CLWall::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Wall" 
      << xmlw::attr("Name") << localName
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Render") << render
      << range
      << xmlw::tag("Norm")
      << vNorm
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << vPosition / Sim->dynamics.units().unitLength()
      << xmlw::endtag("Origin");
}

void 
CLWall::write_povray_info(std::ostream& os) const
{
  const bool asbox = false;
  
  if (render)
    if (asbox)
      os << "object { intersection { object { box { <-0.5, " << -0.5 * Sim->dynamics.units().unitLength() << ", -0.5>, <0.5, " << -0.75 * Sim->dynamics.units().unitLength() << ", 0.5> } Point_At_Trans(<"
	 << vNorm[0] << "," << vNorm[1] << "," << vNorm[2] << ">) translate <"
	 <<  vPosition[0] << "," <<  vPosition[1] << "," <<  vPosition[2] << ">  }\n"
	 << "box { <" 
	 << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
	 << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
	 << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
	 << ">,"
	 << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
	 << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
	 << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
	 << "> }\n"
	 << "} pigment { Col_Glass_Bluish }   }\n";
    else
      {
	Vector pos = vPosition - vNorm * 0.5 * Sim->dynamics.units().unitLength();
	os << "plane { <" << vNorm[0] << "," << vNorm[1] << "," << vNorm[2] 
	   << "> 0 translate <" 
	   <<  pos[0] << "," <<  pos[1] << "," <<  pos[2] 
	   << ">\n"
//Single Colour plane
//	  " texture { pigment { rgb<0.007843137,0.20392156,0.39607843> } } }";
//Tiled effect
" texture {                                                \n\
       tiles {						   \n\
          texture {					   \n\
             pigment {					   \n\
                marble					   \n\
                color_map {				   \n\
                   [0.0 color rgb <0.9, 0.9, 0.9>]	   \n\
                   [0.75 color rgb <0.75, 0.75, 0.75>]	   \n\
                   [1.0 color rgb <0.2, 0.2, 0.2>]	   \n\
                }					   \n\
                turbulence 0.8				   \n\
             }						   \n\
             finish { ambient 1 phong 0.7 }		   \n\
 	    scale 0.5					   \n\
          }						   \n\
          tile2						   \n\
          texture {					   \n\
             pigment {					   \n\
                marble					   \n\
                color_map {				   \n\
                   [0.0 color rgb <0.2, 0.2, 0.2>]         \n\
                   [0.75 color rgb <0.25, 0.25, 0.25>]	   \n\
                   [1.0 color rgb <1, 1, 1>]		   \n\
                }					   \n\
                turbulence 0.8				   \n\
             }						   \n\
             finish { ambient 1 phong 0.7 }			   \n\
             rotate <0, 0, 90>				   \n\
 	    scale 0.25					   \n\
          }						   \n\
       }						   \n\
       scale 0.1					   \n\
    } } 						   \n\
\n";
      }
}

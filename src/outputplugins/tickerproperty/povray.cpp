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

#include "povray.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../base/is_colormap.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/interactions/squarebond.hpp"
#include "../../dynamics/ranges/2RList.hpp"

COPPovray::COPPovray(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPTicker(tmp,"Povray"),
  frameCount(0),
  zoomlevel(1.0)
{
  operator<<(XML);
}

void 
COPPovray::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Zoom"))
	zoomlevel = boost::lexical_cast<Iflt>
	  (XML.getAttribute("Zoom"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPPovray";
    }  
}

void 
COPPovray::ticker()
{
  printImage();
}

void
COPPovray::printImage()
{
  //Dont let this fill up your hard drive!
  if (frameCount > 10000)
    return;

  char *fileName;

  if ( asprintf(&fileName, "%05d", frameCount++) < 0)
    D_throw() << "asprintf error in tinkerXYZ";

  std::ofstream of((std::string("Povray.frame") + fileName + std::string(".pov")).c_str());
  
  free(fileName);

  if (!of.is_open())
    D_throw() << "Could not open povray file for writing";

  //Header of povray file
  of << "#include \"colors.inc\" 	   \n\
#include \"transforms.inc\"                \n\
#declare zoom = "<< zoomlevel << ";	   \n\
#declare cameraloc = <0, zoom, 0>;	   \n\
#declare camera_rot = transform {          \n\
 rotate <20,0,clock*180>                    \n\
};                                          \n\
global_settings { max_trace_level 50 }     \n\
camera {				   \n\
 location cameraloc                        \n\
 transform camera_rot                      \n\
 look_at  <0, 0, 0>			   \n\
}        				   \n\
background { color White }		   \n\
light_source { cameraloc color White transform camera_rot }  \n\
light_source { <2*zoom, 2*zoom, 2*zoom> color White }  \n\
light_source { <-2*zoom, 2*zoom, 2*zoom> color White }  \n\
light_source { <2*zoom, -2*zoom, 2*zoom> color White }  \n\
light_source { <2*zoom, 2*zoom, -2*zoom> color White }  \n\
light_source { <-2*zoom, 2*zoom, -2*zoom> color White }  \n\
light_source { <-2*zoom, -2*zoom, -2*zoom> color White }  \n\
light_source { <-2*zoom, -2*zoom, 2*zoom> color White }  \n\
light_source { <2*zoom, -2*zoom, -2*zoom> color White }  \n\
#include \"glass.inc\"                     \n\
intersection { union {                     \n\
";

  DYNAMO::ColorMap<unsigned int> colmap(0,Sim->Dynamics.getSpecies().size()-1);
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& intPtr, Sim->Dynamics.getInteractions())
    intPtr->write_povray_info(of);

  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->Dynamics.getSpecies())
    spec->getIntPtr()->write_povray_desc
    (colmap.getColor(spec->getID()), spec->getID(), of);
  
  BOOST_FOREACH(const smrtPlugPtr<CLocal>& ptr, Sim->Dynamics.getLocals())
    ptr->write_povray_info(of);


  of << "\n}\nbox { <" 
     << -Sim->aspectRatio[0]/2 - Sim->Dynamics.units().unitLength() 
     << "," << -Sim->aspectRatio[1]/2 - Sim->Dynamics.units().unitLength()  
     << "," << -Sim->aspectRatio[2]/2 - Sim->Dynamics.units().unitLength() 
     << ">,"
     << "<" << Sim->aspectRatio[0]/2 + Sim->Dynamics.units().unitLength()
     << "," << Sim->aspectRatio[1]/2 + Sim->Dynamics.units().unitLength()
     << "," << Sim->aspectRatio[2]/2 + Sim->Dynamics.units().unitLength()
     << "> }\n"
     << "cutaway_textures }\n";
  of.close();
}

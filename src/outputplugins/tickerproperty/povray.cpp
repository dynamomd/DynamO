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

OPPovray::OPPovray(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"Povray"),
  frameCount(0),
  zoomlevel(1.0)
{
  operator<<(XML);
}

void 
OPPovray::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Zoom"))
	zoomlevel = boost::lexical_cast<Iflt>
	  (XML.getAttribute("Zoom"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPPovray";
    }  
}

void 
OPPovray::ticker()
{
  printImage();
}

void
OPPovray::printImage()
{
  //Dont let this fill up your hard drive!
  if (frameCount > 90000)
    return;

  char *fileName;

  const bool boundSystem = false;
  const bool showSky = false;
  const bool fog = false;

  if ( asprintf(&fileName, "%05d", frameCount++) < 0)
    D_throw() << "asprintf error in tinkerXYZ";

  std::ofstream of((std::string("Povray.frame") + fileName + std::string(".pov")).c_str());
  
  free(fileName);

  if (!of.is_open())
    D_throw() << "Could not open povray file for writing";

  //Header of povray file
  of << "\
#include \"colors.inc\" 	           \n\
#include \"transforms.inc\"                \n\
#include \"glass.inc\"                     \n\
global_settings { max_trace_level 20 }     \n\
global_settings { noise_generator 1 }      \n\
global_settings { ambient_light 8 }        \n\
background { rgb<2.0/255.0, 52.0/255.0, 101.0/255.0> }               \n\
// camera ---------------------------------\n\
#declare zoom = "<< zoomlevel << " ;       \n\
#declare Cam0 =                            \n\
   camera {                                \n\
           location  <0.0 , 0 , -zoom>  \n\
           look_at   <0.0 , 0.0 , 0.0>}    \n\
camera{Cam0}                               \n\
// sun ------------------------------------\n\
light_source{<5,1.5,-5> color White}\n\
";

  if (showSky)
    of << "\
// sky ------------------------------------      \n\
global_settings { assumed_gamma 1.0 }            \n\
plane{<0,1,0>,2 hollow                           \n\
      texture{ pigment {color rgb< 0.05,0.25,0.6>} \n\
               finish {ambient 1  diffuse 0}}      \n\
     }                                             \n\
plane{<0,1,0>,1 hollow                           \n\
      texture{pigment { bozo turbulence 0.65                           \n\
                        octaves 6  omega 0.7 lambda 2 		       \n\
                        color_map { [0.0 color rgb <0.95, 0.95, 0.95>] \n\
                                    [0.1 color rgb <0.85, 0.85, 0.85>] \n\
                                    [0.5 color rgbt <1, 1, 1, 1> ]     \n\
                                    [1.0 color rgbt <1, 1, 1, 1> ]   } \n\
                        rotate<10,20,0>				       \n\
                        scale <0.3, 0.4, 0.2>*3 }                     \n\
              finish {ambient 1 diffuse 0}}                            \n\
      }                                                                \n\
";

  if (fog)
    of << "fog { distance 12  color White }                           \n\
fog { distance 2 fog_type 2 fog_alt 0.01 fog_offset -0.1 color White }\n";

  if (boundSystem)
    of << "intersection { union {                     \n";
  

  DYNAMO::ColorMap<unsigned int> colmap(0,Sim->dynamics.getSpecies().size()-1);
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& intPtr, Sim->dynamics.getInteractions())
    intPtr->write_povray_info(of);

  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->dynamics.getSpecies())
    spec->getIntPtr()->write_povray_desc
    (colmap.getColor(spec->getID()), spec->getID(), of);

  BOOST_FOREACH(const smrtPlugPtr<CLocal>& ptr, Sim->dynamics.getLocals())
    ptr->write_povray_info(of);

  //Mirror for the oscillating plate
  //  of << "object { box { <-1, -0.25, 0.01>, <0.5, 0.140977, 0.01> } rotate <45,0,0> translate <0,0,0.30977>  finish { reflection 0.9 ambient 0 diffuse 0 }}";

  if (boundSystem)
    of << "\n}\nbox { <" 
       << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
       << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
       << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
       << ">,"
       << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
       << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
       << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
       << "> }\n"
       << "}\n";
  
  of.close();
}

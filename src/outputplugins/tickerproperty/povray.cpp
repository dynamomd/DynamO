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

COPPovray::COPPovray(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"Povray"),
  frameCount(0)
{}

void 
COPPovray::ticker()
{
  printImage();
}

void
COPPovray::printImage()
{

  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;

  std::string filename = std::string("povray.frame") + boost::lexical_cast<std::string>(frameCount++) + std::string(".pov");
  
  std::ofstream of((std::string("Povray.frame") + boost::lexical_cast<std::string>(frameCount++) + std::string(".pov")).c_str());
  
  if (!of.is_open())
    D_throw() << "Could not open povray file for writing";

  //Header of povray file
  of << "#include \"colors.inc\" 	   \n\
#declare zoom = 1.5;			   \n\
global_settings { max_trace_level 50 }     \n\
camera {				   \n\
 location <0, zoom, 0>			   \n\
 look_at  <0, 0, 0>			   \n\
 rotate <clock*360,clock*180,0>            \n\
}        				   \n\
background { color White }		   \n\
light_source { <0, zoom, 0> color White }  \n\
light_source { <0, -zoom, 0> color White } \n\
light_source { <zoom, 0, 0> color White }  \n\
light_source { <-zoom, 0, 0> color White } \n\
light_source { <0, 0, zoom> color White }  \n\
light_source { <0, 0, -zoom> color White } \n\
";
  DYNAMO::ColorMap<unsigned int> colmap(0,Sim->Dynamics.getSpecies().size()-1);

  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& intPtr, Sim->Dynamics.getInteractions())
    intPtr->write_povray_info(of);

  BOOST_FOREACH(const CSpecies& spec, Sim->Dynamics.getSpecies())
    spec.getIntPtr()->write_povray_desc
    (colmap.getColor(spec.getID()), spec.getID(), of);
  
  of.close();
}

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

#ifdef DYNAMO_visualizer
# include <gtkmm.h>
#endif


#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include "include.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include <boost/tokenizer.hpp>

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include "../liouvillean/CompressionL.hpp"
# include <magnet/HSV.hpp>
# include <coil/RenderObj/Spheres.hpp>
#endif

void 
SpPoint::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
    
    try {
      mass = boost::lexical_cast<double>(XML.getAttribute("Mass"))
	* Sim->dynamics.units().unitMass();
      spName = XML.getAttribute("Name");
      intName = XML.getAttribute("IntName");

      if (XML.isAttributeSet("Color"))
	{
	  typedef boost::tokenizer<boost::char_separator<char> >
	    Tokenizer;
	  
	  boost::char_separator<char> colorSep(",");
	  
	  std::string data(XML.getAttribute("Color"));
	  
	  Tokenizer tokens(data, colorSep);
	  Tokenizer::iterator value_iter = tokens.begin();
	  
	  if (value_iter == tokens.end())
	    throw std::runtime_error("Malformed color in species");
	  _constColor[0] = boost::lexical_cast<int>(*value_iter);
	  
	  if (++value_iter == tokens.end())
	    throw std::runtime_error("Malformed color in species");
	  _constColor[1] = boost::lexical_cast<int>(*value_iter);
	  
	  if (++value_iter == tokens.end())
	    throw std::runtime_error("Malformed color in species");
	  _constColor[2] = boost::lexical_cast<int>(*value_iter);
	  
	  if (++value_iter != tokens.end())
	    throw std::runtime_error("Malformed color in species");
	  
	  _constColor[3] = 255;
	  
	  _colorMode = CONSTANT;
	}

    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in SpPoint";
      }

}

void 
SpPoint::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Mass") 
      << mass / Sim->dynamics.units().unitMass()
      << xml::attr("Name") << spName
      << xml::attr("IntName") << intName
      << xml::attr("Type") << "Point";

  if (_colorMode == CONSTANT)
    {
      std::string colorval = boost::lexical_cast<std::string>(_constColor[0] + 0);
      colorval += ",";
      colorval += boost::lexical_cast<std::string>(_constColor[1] + 0);
      colorval += ",";
      colorval += boost::lexical_cast<std::string>(_constColor[2] + 0);
      XML << xml::attr("Color") << colorval;
    }
  XML << range;
}

void
SpPoint::initialise()
{ 
  if (IntPtr == NULL)
    M_throw() << "SpPoint missing a matching interaction";
}


#ifdef DYNAMO_visualizer

magnet::thread::RefPtr<RenderObj>& 
SpPoint::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      _renderObj = new RTSpheres(range->size(), "Species: " + spName);
      particleData.resize(range->size());
      particleColorData.resize(range->size());
    }

  return _renderObj;
}

void
SpPoint::updateRenderObj(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  //Check if the system is compressing and adjust the radius scaling factor
  float factor = 1;
  if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
    factor = (1 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime);
 
  double diam = getIntPtr()->hardCoreDiam() * factor;
  
  size_t sphID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      Vector pos = Sim->particleList[ID].getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      for (size_t i(0); i < NDIM; ++i)
	particleData[sphID].s[i] = pos[i];
      
      particleData[sphID++].w = diam * 0.5;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getSphereDataBuffer(),
       false, 0, range->size() * sizeof(cl_float4), &particleData[0]);
  }

}

void 
SpPoint::updateColorObj(magnet::CL::CLGLState& CLState) const
{

  switch (_colorMode)
    {
    case IDHSV:
      {
	size_t np = range->size();
	for (size_t sphID(0); sphID < np; ++ sphID)
	  magnet::color::HSVtoRGB(particleColorData[sphID], ((float)(sphID)) / np);
      }
      break;
    case CONSTANT:
      for (size_t id(0); id < range->size(); ++id)
	for (size_t cc(0); cc < 4; ++cc)
	  particleColorData[id].s[cc] = _constColor[cc];
      break;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getColorDataBuffer(),
       false, 0, range->size() * sizeof(cl_uchar4), &particleColorData[0]);
  }
}

void 
SpPoint::showControls(Gtk::ScrolledWindow* win)
{

}
#endif

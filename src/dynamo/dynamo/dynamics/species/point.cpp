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
# include <coil/RenderObj/SphericalParticles.hpp>
# include "../liouvillean/CompressionL.hpp"
#endif

#include "include.hpp"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../units/units.hpp"
#include "../interactions/representations/spherical.hpp"
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

void 
SpPoint::operator<<(const magnet::xml::Node& XML)
{
  range.set_ptr(CRange::getClass(XML,Sim));
  
  try {
    _mass = Sim->_properties.getProperty(XML.getAttribute("Mass"),
					 Property::Units::Mass());
    spName = XML.getAttribute("Name");
    intName = XML.getAttribute("IntName");
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in SpPoint";
    }
}

void 
SpPoint::outputXML(magnet::xml::XmlStream& XML) const
{
  XML << magnet::xml::attr("Mass") 
      << _mass->getName()
      << magnet::xml::attr("Name") << spName
      << magnet::xml::attr("IntName") << intName
      << magnet::xml::attr("Type") << "Point"
      << range;
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
      if (dynamic_cast<const SphericalRepresentation*>(getIntPtr()) == NULL)
	M_throw() << "The interaction " << getIntPtr()->getName() 
		  << " is not able to be drawn using spheres, and yet it is used in the species " << getName()
		  << " as the representative interaction.";

      size_t nsph = dynamic_cast<const SphericalRepresentation&>(*getIntPtr()).spheresPerParticle();

      _renderObj = new RSphericalParticles(nsph * range->size(), "Species: " + spName);
      _coil = new CoilRegister;
    }

  return _renderObj;
}

void
SpPoint::updateRenderData(magnet::GL::Context& context) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  std::vector<cl_float4>& particleData 
    = _renderObj.as<RSphericalParticles>()._particleData;

  ///////////////////////POSITION DATA UPDATE
  //Divide by the maximum box length, to have a natural scale for the visualizer
  const double lengthRescale = 1 / Sim->primaryCellSize.maxElement();
  
  //Check if the system is compressing and adjust the radius scaling factor
  float rfactor = lengthRescale;
  if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
    rfactor *= (1 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime);
  
  const SphericalRepresentation& data
    = dynamic_cast<const SphericalRepresentation&>(*getIntPtr());

  size_t nsph = data.spheresPerParticle();

  size_t sphID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      for (size_t s(0); s < nsph; ++s)
	{
	  Vector pos = data.getPosition(ID, s);
	  
	  for (size_t i(0); i < NDIM; ++i)
	    particleData[s * range->size() + sphID].s[i] = pos[i] * lengthRescale;
	  
	  particleData[s * range->size() + sphID].w = 0.5 * rfactor * data.getDiameter(ID, s);
	}
      ++sphID;
    }

  _renderObj.as<RSphericalParticles>().recolor(context);
  
  _coil->getInstance().getTaskQueue().queueTask(magnet::function::Task::makeTask
						(&RSphericalParticles::sendRenderData,
						 &(_renderObj.as<RSphericalParticles>()), 
						 context));
}
#endif

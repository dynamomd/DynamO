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

#include "include.hpp"
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

#ifdef DYNAMO_visualizer
#include "../interactions/representations/spherical.hpp"
#include "../liouvillean/CompressionL.hpp"
#endif

Species* 
Species::getClass(const magnet::xml::Node& XML, dynamo::SimData* tmp, size_t nID)
{
  if (!std::strcmp(XML.getAttribute("Type"), "Point"))
    return new SpPoint(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "SphericalTop"))
    return new SpSphericalTop(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "Lines"))
    return new SpLines(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "Dumbbells"))
    return new SpDumbbells(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "FixedCollider"))
    return new SpFixedCollider(XML, tmp, nID);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of species encountered";
}


#ifdef DYNAMO_visualizer

std::tr1::shared_ptr<coil::DataSet>
Species::createDataSet() const
{
  if (dynamic_cast<const SphericalRepresentation*>(getIntPtr()) == NULL)
    M_throw() << "The interaction " << getIntPtr()->getName() 
	      << " is not able to be drawn using spheres, and yet it is used in the species " << getName()
	      << " as the representative interaction.";
  
  size_t nsph = dynamic_cast<const SphericalRepresentation&>(*getIntPtr()).spheresPerParticle();

  _renderData.reset(new coil::DataSet("Species: " + spName, nsph * range->size()));
  return _renderData;
}

void
Species::initDataSet() const
{
  _renderData->addAttribute("Positions", coil::Attribute::COORDINATE, 3);
  _renderData->addAttribute("Mass", coil::Attribute::EXTENSIVE, 1);
  { 
    size_t nsph = dynamic_cast<const SphericalRepresentation&>(*getIntPtr()).spheresPerParticle();    
    std::vector<GLfloat>& mass = (*_renderData)["Mass"].getData();
    size_t sphID(0);
    BOOST_FOREACH(unsigned long ID, *range)
      {
	for (size_t s(0); s < nsph; ++s)
	  mass[nsph * sphID + s] = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
	++sphID;
      }
    (*_renderData)["Mass"].flagNewData();
  }
  _renderData->addAttribute("Radii", coil::Attribute::INTENSIVE, 1);
}

void
Species::updateRenderData() const
{
  if (!_renderData)
    M_throw() << "Updating before the render object has been fetched";
  
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

  std::vector<GLfloat>& posdata = (*_renderData)["Positions"].getData();
  std::vector<GLfloat>& radii = (*_renderData)["Radii"].getData();

  size_t sphID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      for (size_t s(0); s < nsph; ++s)
	{
	  Vector pos = data.getPosition(ID, s);
	  
	  for (size_t i(0); i < NDIM; ++i)
	    posdata[3 * (nsph * sphID + s) + i] = pos[i] * lengthRescale;
	  
	  radii[nsph * sphID + s] = 0.5 * rfactor * data.getDiameter(ID, s);
	}

      ++sphID;
    }

  (*_renderData)["Positions"].flagNewData();
  (*_renderData)["Radii"].flagNewData();
}
#endif

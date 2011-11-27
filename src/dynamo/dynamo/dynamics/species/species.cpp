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

#include <dynamo/dynamics/species/include.hpp>

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/DataSet.hpp>
# include <dynamo/dynamics/interactions/representations/spherical.hpp>
# include <dynamo/dynamics/liouvillean/CompressionL.hpp>
#endif

#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  Species::~Species() {}

  shared_ptr<Species>
  Species::getClass(const magnet::xml::Node& XML, dynamo::SimData* tmp, size_t nID)
  {
    if (!std::strcmp(XML.getAttribute("Type"), "Point"))
      return shared_ptr<Species>(new SpPoint(XML, tmp, nID));
    else if (!std::strcmp(XML.getAttribute("Type"), "SphericalTop"))
      return shared_ptr<Species>(new SpSphericalTop(XML, tmp, nID));
    else if (!std::strcmp(XML.getAttribute("Type"), "Lines"))
      return shared_ptr<Species>(new SpLines(XML, tmp, nID));
    else if (!std::strcmp(XML.getAttribute("Type"), "Dumbbells"))
      return shared_ptr<Species>(new SpDumbbells(XML, tmp, nID));
    else if (!std::strcmp(XML.getAttribute("Type"), "FixedCollider"))
      return shared_ptr<Species>(new SpFixedCollider(XML, tmp, nID));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of species encountered";
  }


#ifdef DYNAMO_visualizer

  shared_ptr<coil::DataSet>
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
    _renderData->addAttribute("Positions", coil::Attribute::COORDINATE | coil::Attribute::DEFAULT_GLYPH_POSITION, 3);
    _renderData->addAttribute("Velocity", coil::Attribute::INTENSIVE, 3);
    _renderData->addAttribute("Radii", coil::Attribute::INTENSIVE | coil::Attribute::DEFAULT_GLYPH_SCALING, 1);
    _renderData->addAttribute("Mass", coil::Attribute::EXTENSIVE, 1);

    { 
      size_t nsph = dynamic_cast<const SphericalRepresentation&>(*getIntPtr()).spheresPerParticle();    
      std::vector<GLfloat>& mass = (*_renderData)["Mass"];
      size_t sphID(0);
      BOOST_FOREACH(unsigned long ID, *range)
	{
	  for (size_t s(0); s < nsph; ++s)
	    mass[nsph * sphID + s] = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
	  ++sphID;
	}
      (*_renderData)["Mass"].flagNewData();
    }

    _renderData->addAttribute("ID", coil::Attribute::INTENSIVE, 1);
    { 
      size_t nsph = dynamic_cast<const SphericalRepresentation&>(*getIntPtr()).spheresPerParticle();    
      std::vector<GLfloat>& mass = (*_renderData)["ID"];
      size_t sphID(0);
      BOOST_FOREACH(unsigned long ID, *range)
	{
	  for (size_t s(0); s < nsph; ++s)
	    mass[nsph * sphID + s] = ID;
	  ++sphID;
	}
      (*_renderData)["ID"].flagNewData();
    }
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

    std::vector<GLfloat>& posdata = (*_renderData)["Positions"];
    std::vector<GLfloat>& veldata = (*_renderData)["Velocity"];
    std::vector<GLfloat>& radii = (*_renderData)["Radii"];

    size_t sphID(0);
    BOOST_FOREACH(unsigned long ID, *range)
      {
	Vector vel = Sim->particleList[ID].getVelocity();
	for (size_t s(0); s < nsph; ++s)
	  {
	    Vector pos = data.getPosition(ID, s);
	  
	    for (size_t i(0); i < NDIM; ++i)
	      posdata[3 * (nsph * sphID + s) + i] = pos[i] * lengthRescale;

	    for (size_t i(0); i < NDIM; ++i)
	      veldata[3 * (nsph * sphID + s) + i] = vel[i] * lengthRescale;
	  
	    radii[nsph * sphID + s] = 0.5 * rfactor * data.getDiameter(ID, s);
	  }
	++sphID;
      }

    (*_renderData)["Positions"].flagNewData();
    (*_renderData)["Velocity"].flagNewData();
    (*_renderData)["Radii"].flagNewData();
  }
#endif
}

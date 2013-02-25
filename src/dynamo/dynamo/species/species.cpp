/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/species/include.hpp>

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/DataSet.hpp>
# include <dynamo/interactions/glyphrepresentation.hpp>
# include <dynamo/dynamics/compression.hpp>
# include <dynamo/schedulers/scheduler.hpp>
# include <dynamo/BC/LEBC.hpp>
#endif

#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  Species::~Species() {}

  shared_ptr<Species>
  Species::getClass(const magnet::xml::Node& XML, dynamo::Simulation* tmp, size_t nID)
  {
    if (!XML.getAttribute("Type").getValue().compare("Point"))
      return shared_ptr<Species>(new SpPoint(XML, tmp, nID));
    else if (!XML.getAttribute("Type").getValue().compare("SphericalTop")
	     || !XML.getAttribute("Type").getValue().compare("Lines"))
      return shared_ptr<Species>(new SpSphericalTop(XML, tmp, nID));
    else if (!XML.getAttribute("Type").getValue().compare("FixedCollider"))
      return shared_ptr<Species>(new SpFixedCollider(XML, tmp, nID));
    else 
      M_throw() << XML.getAttribute("Type").getValue()
		<< ", Unknown type of species encountered";
  }


#ifdef DYNAMO_visualizer

  shared_ptr<coil::DataSet>
  Species::createDataSet() const
  {
    if (dynamic_cast<const GlyphRepresentation*>(getIntPtr()) == NULL)
      M_throw() << "The interaction " << getIntPtr()->getName() 
		<< " is not able to be drawn by the visualiser, and yet it is used in the species " << getName()
		<< " as the representative interaction.";
  
    const GlyphRepresentation& representation = dynamic_cast<const GlyphRepresentation&>(*getIntPtr());
    size_t nsph = representation.glyphsPerParticle();
    
    _renderData.reset(new coil::DataSet("Species: " + spName, nsph * range->size(), representation.getDefaultGlyphType()));
    return _renderData;
  }

  void
  Species::initDataSet() const
  {
    _renderData->addAttribute("Position", coil::Attribute::COORDINATE | coil::Attribute::DEFAULT_GLYPH_POSITION, 3);
    _renderData->addAttribute("Velocity", coil::Attribute::INTENSIVE, 3);
    _renderData->addAttribute("Size", coil::Attribute::INTENSIVE | coil::Attribute::DEFAULT_GLYPH_SCALING, 3);
    _renderData->addAttribute("Mass", coil::Attribute::EXTENSIVE, 1);
    _renderData->addAttribute("Event Count", coil::Attribute::EXTENSIVE, 1);

    _renderData->setPeriodicVectors(Vector(Sim->primaryCellSize[0], 0, 0),
				    Vector(0, Sim->primaryCellSize[1], 0),
				    Vector(0, 0, Sim->primaryCellSize[2]));

    if (Sim->dynamics->hasOrientationData())
      {
	_renderData->addAttribute("Orientation", coil::Attribute::EXTENSIVE | coil::Attribute::DEFAULT_GLYPH_ORIENTATION, 3);
	_renderData->addAttribute("Angular Velocity", coil::Attribute::EXTENSIVE, 3);
      }

    { 
      size_t nsph = dynamic_cast<const GlyphRepresentation&>(*getIntPtr()).glyphsPerParticle();    
      std::vector<GLfloat>& mass = (*_renderData)["Mass"];
      size_t sphID(0);
      BOOST_FOREACH(unsigned long ID, *range)
	{
	  for (size_t s(0); s < nsph; ++s)
	    mass[nsph * sphID + s] = Sim->species[Sim->particles[ID]]->getMass(ID) / Sim->units.unitMass();
	  ++sphID;
	}
      (*_renderData)["Mass"].flagNewData();
    }

    _renderData->addAttribute("ID", coil::Attribute::INTENSIVE, 1);
    { 
      size_t nsph = dynamic_cast<const GlyphRepresentation&>(*getIntPtr()).glyphsPerParticle();    
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
    
    _renderData->getContext()->queueTask(magnet::function::Task::makeTask(&coil::DataSet::addGlyphs, _renderData.get()));
  }

  void
  Species::updateRenderData() const
  {
    if (!_renderData)
      M_throw() << "Updating before the render object has been fetched";
    
    shared_ptr<BCLeesEdwards> BC = std::tr1::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs);
    if (BC)
      {
	_renderData->setPeriodicVectors(Vector(Sim->primaryCellSize[0], 0, 0),
					Vector(BC->getBoundaryDisplacement(), Sim->primaryCellSize[1], 0),
					Vector(0, 0, Sim->primaryCellSize[2]));
      }


    ///////////////////////POSITION DATA UPDATE
    //Check if the system is compressing and adjust the radius scaling factor
    float rfactor = 1.0;
    if (std::tr1::dynamic_pointer_cast<DynCompression>(Sim->dynamics))
      rfactor *= (1 + static_cast<const DynCompression&>(*Sim->dynamics).getGrowthRate() * Sim->systemTime);
  
    const GlyphRepresentation& data
      = dynamic_cast<const GlyphRepresentation&>(*getIntPtr());

    const size_t nsph = data.glyphsPerParticle();
    std::vector<GLfloat>& posdata = (*_renderData)["Position"];
    std::vector<GLfloat>& veldata = (*_renderData)["Velocity"];
    std::vector<GLfloat>& sizes = (*_renderData)["Size"];
    std::vector<GLfloat>& eventCounts = (*_renderData)["Event Count"];
    const std::vector<size_t>& simEventCounts = Sim->ptrScheduler->getEventCounts();
    
    size_t glyphID(0);
    BOOST_FOREACH(unsigned long ID, *range)
      {
	Vector vel = Sim->particles[ID].getVelocity() / Sim->units.unitVelocity();
	for (size_t s(0); s < nsph; ++s)
	  {
	    Vector pos = data.getGlyphPosition(ID, s) / Sim->units.unitLength();
	  
	    for (size_t i(0); i < NDIM; ++i)
	      posdata[3 * (nsph * glyphID + s) + i] = pos[i];

	    for (size_t i(0); i < NDIM; ++i)
	      veldata[3 * (nsph * glyphID + s) + i] = vel[i];

	    Vector psize = rfactor * data.getGlyphSize(ID, s) / Sim->units.unitLength();
	    for (size_t i(0); i < NDIM; ++i)
	      sizes[3 * (nsph * glyphID + s) + i] = psize[i];
	  }

	eventCounts[glyphID] = 0;
	if (!simEventCounts.empty())
	  eventCounts[glyphID] = simEventCounts[ID];
	++glyphID;
      }

    if (Sim->dynamics->hasOrientationData())
      {
	std::vector<GLfloat>& orientationdata = (*_renderData)["Orientation"];
	std::vector<GLfloat>& angularvdata = (*_renderData)["Angular Velocity"];
	size_t glyphID(0);
	const std::vector<Dynamics::rotData>& data = Sim->dynamics->getCompleteRotData();
	BOOST_FOREACH(unsigned long ID, *range)
	  {
	    for (size_t s(0); s < nsph; ++s)
	      {
		for (size_t i(0); i < NDIM; ++i)
		  {
		    orientationdata[3 * (nsph * glyphID + s) + i] = data[ID].orientation[i];
		    angularvdata[3 * (nsph * glyphID + s) + i] = data[ID].angularVelocity[i] * Sim->units.unitTime();
		  }
	      }
	    ++glyphID;
	  }
	(*_renderData)["Angular Velocity"].flagNewData();
	(*_renderData)["Orientation"].flagNewData();
      }

    (*_renderData)["Position"].flagNewData();
    (*_renderData)["Velocity"].flagNewData();
    (*_renderData)["Size"].flagNewData();
    (*_renderData)["Event Count"].flagNewData();
  }
#endif
}

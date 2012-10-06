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

#include <dynamo/globals/cells.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/ranges/1RAll.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/ranges/1RList.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstdio>
#include <set>
#include <algorithm>

namespace dynamo {
  GCells::GCells(dynamo::Simulation* nSim, const std::string& name, size_t overlink):
    GNeighbourList(nSim, "CellNeighbourList"),
    cellDimension(1,1,1),
    _oversizeCells(1.0),
    NCells(0),
    overlink(overlink)
  {
    globName = name;
    dout << "Cells Loaded" << std::endl;
  }

  GCells::GCells(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    GNeighbourList(ptrSim, "CellNeighbourList"),
    cellDimension(1,1,1),
    _oversizeCells(1.0),
    NCells(0),
    overlink(1)
  {
    operator<<(XML);

    dout << "Cells Loaded" << std::endl;
  }

  void 
  GCells::operator<<(const magnet::xml::Node& XML)
  {
    try {
      if (XML.hasAttribute("OverLink"))
	overlink = XML.getAttribute("OverLink").as<size_t>();

      if (XML.hasAttribute("NeighbourhoodRange"))
	_maxInteractionRange = XML.getAttribute("NeighbourhoodRange").as<double>()
	  * Sim->units.unitLength();

      if (XML.hasAttribute("Oversize"))
	_oversizeCells = XML.getAttribute("Oversize").as<double>();

      if (_oversizeCells < 1.0)
	M_throw() << "You must specify an Oversize greater than 1.0, otherwise your cells are too small!";
    
      globName = XML.getAttribute("Name");

      if (XML.hasAttribute("Range"))
	range = shared_ptr<Range>(Range::getClass(XML, Sim));
    }
    catch(...)
      {
	M_throw() << "Error loading GCells";
      }
  }

  GlobalEvent 
  GCells::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    //This 
    //Sim->dynamics->updateParticle(part);
    //is not required as we compensate for the delay using 
    //Sim->dynamics->getParticleDelay(part)
  
    return GlobalEvent(part,
		       Sim->dynamics->
		       getSquareCellCollision2
		       (part, 
			calcPosition(partCellData[part.getID()], part), 
			cellDimension)
		       -Sim->dynamics->getParticleDelay(part),
		       CELL, *this);

  }

  void
  GCells::runEvent(Particle& part, const double) const
  {
    //Despite the system not being streamed this must be done.  This is
    //because the scheduler and all interactions, locals and systems
    //expect the particle to be up to date.
    Sim->dynamics->updateParticle(part);

    boost::unordered_map<size_t, size_t>::iterator it = partCellData.find(part.getID());

    const size_t oldCell(it->second);

    size_t endCell;

    //Determine the cell transition direction, its saved
    int cellDirectionInt(Sim->dynamics->
			 getSquareCellCollision3
			 (part, calcPosition(oldCell, part), cellDimension));
  
    size_t cellDirection = abs(cellDirectionInt) - 1;

    //The coordinates of the new center cell in the neighbourhood of the
    //particle
    magnet::math::MortonNumber<3> newNBCell(oldCell);

    {
      //The position of the cell the particle will end up in
      magnet::math::MortonNumber<3> dendCell(newNBCell);
    
      if (cellDirectionInt > 0)
	{
	  dendCell[cellDirection] = (dendCell[cellDirection].getRealValue() + 1) % cellCount[cellDirection];
	  newNBCell[cellDirection] = (dendCell[cellDirection].getRealValue() + overlink) % cellCount[cellDirection];
	}
      else
	{
	  //We use the trick of adding cellCount to convert the
	  //subtraction to an addition, to prevent errors in the modulus
	  //of underflowing unsigned integers.
	  dendCell[cellDirection] = (dendCell[cellDirection].getRealValue() 
				     + cellCount[cellDirection] - 1) % cellCount[cellDirection];
	  newNBCell[cellDirection] = (dendCell[cellDirection].getRealValue() 
				      + cellCount[cellDirection] - overlink) % cellCount[cellDirection];
	}
      endCell = dendCell.getMortonNum();
    }

    removeFromCellwIt(part.getID(), it);
    addToCell(part.getID(), endCell);

    //Get rid of the virtual event we're running, an updated event is
    //pushed after all other events are added
    Sim->ptrScheduler->popNextEvent();


    //Particle has just arrived into a new cell warn the scheduler about
    //its new neighbours so it can add them to the heap
    //Holds the displacement in each dimension, the unit is cells!

    //These are the two dimensions to walk in
    size_t dim1 = (cellDirection + 1) % 3,
      dim2 = (cellDirection + 2) % 3;

    newNBCell[dim1] += cellCount[dim1] - overlink;
    newNBCell[dim2] += cellCount[dim2] - overlink;

    const magnet::math::DilatedInteger<3> saved_coord(newNBCell[dim1]);

    //We now have the lowest cell coord, or corner of the cells to update
    for (size_t iDim(0); iDim < 2 * overlink + 1; ++iDim)
      {
	newNBCell[dim2] %= cellCount[dim2];

	for (size_t jDim(0); jDim < 2 * overlink + 1; ++jDim)
	  {
	    newNBCell[dim1] %= cellCount[dim1];
	    
	    BOOST_FOREACH(const size_t& next, list[newNBCell.getMortonNum()])
	      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	        nbs.second(part, next);
	  
	    ++newNBCell[dim1];
	  }

	newNBCell[dim1] = saved_coord; 
	++newNBCell[dim2];
      }

    //Tell about the new locals
    BOOST_FOREACH(const size_t& lID, cells[endCell])
      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewLocalNotify)
        nbs.second(part, lID);
  
    //Push the next virtual event, this is the reason the scheduler
    //doesn't need a second callback
    Sim->ptrScheduler->pushEvent(part, getEvent(part));
    Sim->ptrScheduler->sort(part);

    BOOST_FOREACH(const nbHoodSlot& nbs, sigCellChangeNotify)
      nbs.second(part, oldCell);
  
    //Debug section
#ifdef DYNAMO_WallCollDebug
    {
      magnet::math::MortonNumber<3> newNBCellv(oldCell);
      magnet::math::MortonNumber<3> endCellv(endCell);
    
      dout << "CellEvent: t=" 
	   << Sim->systemTime / Sim->units.unitTime()
	   << " ID "
	   << part.getID()
	   << "  from <" 
	   << newNBCellv[0].getRealValue() << "," << newNBCellv[1].getRealValue() 
	   << "," << newNBCellv[2].getRealValue()
	   << "> to <" 
	   << endCellv[0].getRealValue() << "," << endCellv[1].getRealValue() 
	   << "," << endCellv[2].getRealValue() << ">"
	   << std::endl;
    }
#endif
  }

  void 
  GCells::initialise(size_t nID)
  {
    ID=nID;
    typedef void (GCells::*CallBackType)(size_t);

    _particleAdded = Sim->particle_added_signal()
      .connect(boost::bind(CallBackType(&GCells::addToCell), this, _1));
    _particleRemoved = Sim->particle_removed_signal()
      .connect(boost::bind(CallBackType(&GCells::removeFromCell), this, _1));

    reinitialise();

    dout << "Neighbourlist contains " << partCellData.size() 
	 << " particle entries";
  }

  void
  GCells::reinitialise()
  {
    GNeighbourList::reinitialise();
      
    dout << "Reinitialising on collision " << Sim->eventCount << std::endl;

    //Create the cells
    addCells((_maxInteractionRange 
	      * (1.0 + 10 * std::numeric_limits<double>::epsilon()))
	     * _oversizeCells / overlink);

    addLocalEvents();

    BOOST_FOREACH(const initSlot& nbs, sigReInitNotify)
      nbs.second();

    if (isUsedInScheduler)
      Sim->ptrScheduler->initialise();
  }

  void
  GCells::outputXML(magnet::xml::XmlStream& XML, const std::string& type) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << type
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("NeighbourhoodRange") 
	<< _maxInteractionRange / Sim->units.unitLength();

    if (overlink > 1)   XML << magnet::xml::attr("OverLink") << overlink;
    if (_oversizeCells != 1.0) XML << magnet::xml::attr("Oversize") << _oversizeCells;
  
    XML << *range
	<< magnet::xml::endtag("Global");
  }

  void
  GCells::outputXML(magnet::xml::XmlStream& XML) const
  { outputXML(XML, "Cells"); }

  void
  GCells::addCells(double maxdiam)
  {
    cells.clear();
    list.clear();
    partCellData.clear();
    NCells = 1;

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      {
	cellCount[iDim] = int(Sim->primaryCellSize[iDim] 
			      / (maxdiam * (1.0 + 10 * std::numeric_limits<double>::epsilon())));
      
	if (cellCount[iDim] < 2 * overlink + 1)
	  cellCount[iDim] = 2 * overlink + 1;
	
	NCells *= cellCount[iDim];
      
	dilatedCellMax[iDim] = cellCount[iDim] - 1;
	cellLatticeWidth[iDim] = Sim->primaryCellSize[iDim] / cellCount[iDim];
	cellDimension[iDim] = cellLatticeWidth[iDim] + (cellLatticeWidth[iDim] - maxdiam) 
	  * lambda;
	cellOffset[iDim] = -(cellLatticeWidth[iDim] - maxdiam) * lambda * 0.5;
      }

    dout << "Cells <x,y,z> " << cellCount[0] << ","
	 << cellCount[1] << "," << cellCount[2]
	 << "\nCell Offset "
	 << cellOffset[0] / Sim->units.unitLength() << ","
	 << cellOffset[1] / Sim->units.unitLength() << ","
	 << cellOffset[2] / Sim->units.unitLength()
	 << "\nCells Dimension " 
	 << cellDimension[0] / Sim->units.unitLength()
	 << ","
	 << cellDimension[1] / Sim->units.unitLength()
	 << "," 
	 << cellDimension[2] / Sim->units.unitLength()
	 << "\nLattice spacing " 
	 << cellLatticeWidth[0] / Sim->units.unitLength()
	 << ","
	 << cellLatticeWidth[1] / Sim->units.unitLength()
	 << "," 
	 << cellLatticeWidth[2] / Sim->units.unitLength()
	 << "\nRequested supported length " << maxdiam / Sim->units.unitLength()
	 << "\nSupported length           " << getMaxSupportedInteractionLength() / Sim->units.unitLength()
	 << std::endl;

    if (getMaxSupportedInteractionLength() < maxdiam)
      M_throw() << "The system size is too small to support the range of interactions specified (i.e. the system is smaller than the interaction diameter of one particle).";

    //Find the required size of the morton array
    magnet::math::MortonNumber<3> coords(cellCount[0], cellCount[1], cellCount[2]);
    size_t sizeReq = coords.getMortonNum();

    cells.resize(sizeReq); //Empty Cells created!
    list.resize(sizeReq); //Empty Cells created!

    dout << "Vector Size <N>  " << sizeReq << std::endl;
  
    //Add the particles section
    //Required so particles find the right owning cell
    Sim->dynamics->updateAllParticles();
  
    ////Add all the particles 
    BOOST_FOREACH(const size_t& id, *range)
      {
	Particle& p = Sim->particles[id];
	Sim->dynamics->updateParticle(p); 
	addToCell(id);
#ifdef DYNAMO_WallCollDebug
	boost::unordered_map<size_t, size_t>::iterator it = partCellData.find(id);
	magnet::math::MortonNumber<3> currentCell(it->second);
	
	dout << "Added particle ID=" << id << " to cell <"
	     << currentCell[0].getRealValue() 
	     << "," << currentCell[1].getRealValue()
	     << "," << currentCell[2].getRealValue()
	     << ">" << std::endl;
#endif
      }

    dout << "Cell loading " << float(partCellData.size()) / NCells 
	 << std::endl;
  }

  void 
  GCells::addLocalEvents()
  {
    for (size_t iDim = 0; iDim < cellCount[0]; ++iDim)
      for (size_t jDim = 0; jDim < cellCount[1]; ++jDim)
	for (size_t kDim = 0; kDim < cellCount[2]; ++kDim)
	  {
	    magnet::math::MortonNumber<3> coords(iDim, jDim, kDim);
	    size_t id = coords.getMortonNum();
	    cells[id].clear();
	    Vector pos = calcPosition(coords);
	  
	    //We make the box slightly larger to ensure objects on the boundary are included
	    BOOST_FOREACH(const shared_ptr<Local>& local, Sim->locals)
	      if (local->isInCell(pos - 0.0001 * cellDimension, 1.0002 * cellDimension))
		cells[id].push_back(local->getID());
	  }
  }

  magnet::math::MortonNumber<3>
  GCells::getCellID(Vector pos) const
  {
    Sim->BCs->applyBC(pos);

    magnet::math::MortonNumber<3> retval;

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      {
	int coord = std::floor((pos[iDim] + 0.5 * Sim->primaryCellSize[iDim] - cellOffset[iDim])
			       / cellLatticeWidth[iDim]);
	coord %= cellCount[iDim];
	if (coord < 0) coord += cellCount[iDim];
	retval[iDim] = coord;
      }

    return retval;
  }

  RList
  GCells::getParticleNeighbours(const Particle& part) const
  {
    const magnet::math::MortonNumber<3> particle_cell_coords(partCellData[part.getID()]);

    magnet::math::MortonNumber<3> zero_coords;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      zero_coords[iDim] = (particle_cell_coords[iDim].getRealValue() + cellCount[iDim] - overlink)
	% cellCount[iDim];
    
    RList retval;
    //This initial reserve greatly speeds up the later inserts
    retval.getContainer().reserve(32);

    magnet::math::MortonNumber<3> coords(zero_coords);
    for (size_t x(0); x < 2 * overlink + 1; ++x)
      {
	coords[0] = (zero_coords[0].getRealValue() + x) % cellCount[0];
	for (size_t y(0); y < 2 * overlink + 1; ++y)
	  {
	    coords[1] = (zero_coords[1].getRealValue() + y) % cellCount[1];
	    for (size_t z(0); z < 2 * overlink + 1; ++z)
	      {
		coords[2] = (zero_coords[2].getRealValue() + z) % cellCount[2];

		const std::vector<size_t>&  nlist = list[coords.getMortonNum()];
		retval.getContainer().insert(retval.getContainer().end(), nlist.begin(), nlist.end());
	      }
	  }
      }

    return retval;
  }

  RList
  GCells::getParticleLocals(const Particle& part) const
  {
    return RList(cells[partCellData[part.getID()]]);
  }

  void 
  GCells::getParticleNeighbourhood(const Particle& part,
				   const nbHoodFunc& func) const
  {
    const magnet::math::MortonNumber<3> particle_cell_coords(partCellData[part.getID()]);

    magnet::math::MortonNumber<3> zero_coords;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      zero_coords[iDim] = (particle_cell_coords[iDim].getRealValue() + cellCount[iDim] - overlink)
	% cellCount[iDim];

    magnet::math::MortonNumber<3> coords(zero_coords);
    for (size_t x(0); x < 2 * overlink + 1; ++x)
      {
	coords[0] = (zero_coords[0].getRealValue() + x) % cellCount[0];
	for (size_t y(0); y < 2 * overlink + 1; ++y)
	  {
	    coords[1] = (zero_coords[1].getRealValue() + y) % cellCount[1];
	    for (size_t z(0); z < 2 * overlink + 1; ++z)
	      {
		coords[2] = (zero_coords[2].getRealValue() + z) % cellCount[2];

		BOOST_FOREACH(const size_t& next, list[coords.getMortonNum()])
		  func(part, next);
	      }
	  }
      }
  }

  void
  GCells::getParticleNeighbourhood(const Vector& vec,
				   const nbHoodFunc2& func) const
  {
    const magnet::math::MortonNumber<3> particle_cell_coords(getCellID(vec));
    
    magnet::math::MortonNumber<3> zero_coords;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      zero_coords[iDim] = (particle_cell_coords[iDim].getRealValue() + cellCount[iDim] - overlink)
	% cellCount[iDim];

    magnet::math::MortonNumber<3> coords(zero_coords);
    for (size_t x(0); x < 2 * overlink + 1; ++x)
      {
	coords[0] = (zero_coords[0].getRealValue() + x) % cellCount[0];
	for (size_t y(0); y < 2 * overlink + 1; ++y)
	  {
	    coords[1] = (zero_coords[1].getRealValue() + y) % cellCount[1];
	    for (size_t z(0); z < 2 * overlink + 1; ++z)
	      {
		coords[2] = (zero_coords[2].getRealValue() + z) % cellCount[2];

		BOOST_FOREACH(const size_t& next, list[coords.getMortonNum()])
		  func(next);
	      }
	  }
      }
  }

  void 
  GCells::getLocalNeighbourhood(const Particle& part, 
				const nbHoodFunc& func) const
  {
    BOOST_FOREACH(const size_t& id, cells[partCellData[part.getID()]])
      func(part, id);
  }

  double 
  GCells::getMaxSupportedInteractionLength() const
  {
    double retval(HUGE_VAL);

    for (size_t i = 0; i < NDIM; ++i)
      {
	double supported_length = cellLatticeWidth[i]
	  + lambda * (cellLatticeWidth[i] - cellDimension[i]);

	//Test if, in this dimension, one neighbourhood of cells spans
	//the system. If so, the maximum interaction supported is the
	//system width.
	if (cellCount[i] == 2 * overlink + 1)
	  supported_length = Sim->primaryCellSize[i];

	retval = std::min(retval, supported_length);
      }

    return retval;
  }

  Vector 
  GCells::calcPosition(const magnet::math::MortonNumber<3>& coords, const Particle& part) const
  {
    //We always return the cell that is periodically nearest to the particle
    Vector primaryCell;
  
    for (size_t i(0); i < NDIM; ++i)
      primaryCell[i] = coords[i].getRealValue() * cellLatticeWidth[i]
	- 0.5 * Sim->primaryCellSize[i] + cellOffset[i];

  
    Vector imageCell;
  
    for (size_t i = 0; i < NDIM; ++i)
      imageCell[i] = primaryCell[i]
	- Sim->primaryCellSize[i] * lrint((primaryCell[i] - part.getPosition()[i]) 
					  / Sim->primaryCellSize[i]);

    return imageCell;
  }

  Vector 
  GCells::calcPosition(const magnet::math::MortonNumber<3>& coords) const
  {
    Vector primaryCell;
  
    for (size_t i(0); i < NDIM; ++i)
      primaryCell[i] = coords[i].getRealValue() * cellLatticeWidth[i] 
	- 0.5 * Sim->primaryCellSize[i] + cellOffset[i];
  
    return primaryCell;
  }
}

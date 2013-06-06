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
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/ranges/IDRangeList.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstdio>
#include <set>
#include <algorithm>

namespace {
  const bool verbose = false;
}

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
    
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
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

    if (verbose)
      {
	Vector cellPos = calcPosition(partCellData[part.getID()], part);
	Vector relpos = part.getPosition() - cellPos;
	Sim->BCs->applyBC(relpos);
	derr 
	  << "Calculating event for particle " << part.getID() << " in Cell " << magnet::math::MortonNumber<3>(partCellData[part.getID()]).toString()
	  << "\nParticle pos = " << part.getPosition().toString()
	  << "\nCell pos = " << cellPos.toString()
	  << "\nRelpos = " << relpos.toString()
	  << "\nCell size = " << cellDimension.toString()
	  << "\nTime = " << Sim->dynamics->getSquareCellCollision2(part, calcPosition(partCellData[part.getID()], part),
								   cellDimension) - Sim->dynamics->getParticleDelay(part)
	  << "\nDelay = " << Sim->dynamics->getParticleDelay(part)
	  << std::endl;
      }
  
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

    std::unordered_map<size_t, size_t>::iterator it = partCellData.find(part.getID());

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
    //pushed after the callbacks are complete (the callbacks may also
    //add events so this must be done first).
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
	    
	    for (const size_t& next : list[newNBCell.getMortonNum()])
	      _sigNewNeighbour(part, next);
	  
	    ++newNBCell[dim1];
	  }

	newNBCell[dim1] = saved_coord; 
	++newNBCell[dim2];
      }
  
    //Push the next virtual event, this is the reason the scheduler
    //doesn't need a second callback
    Sim->ptrScheduler->pushEvent(part, getEvent(part));
    Sim->ptrScheduler->sort(part);

    _sigCellChange(part, oldCell);
  
    if (verbose)
      {
	magnet::math::MortonNumber<3> newNBCellv(oldCell);
	magnet::math::MortonNumber<3> endCellv(endCell);
    
	derr << "CellEvent: t=" 
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
  }

  void 
  GCells::initialise(size_t nID)
  {
    ID=nID;

    reinitialise();

    dout << "Neighbourlist contains " << partCellData.size() 
	 << " particle entries"
	 << std::endl;
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

    _sigReInitialise();

    if (isUsedInScheduler)
      Sim->ptrScheduler->initialise();
  }

  void
  GCells::outputXML(magnet::xml::XmlStream& XML) const
  { 
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "Cells"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("NeighbourhoodRange") 
	<< _maxInteractionRange / Sim->units.unitLength();
    
    if (overlink > 1)   XML << magnet::xml::attr("OverLink") << overlink;
    if (_oversizeCells != 1.0) XML << magnet::xml::attr("Oversize") << _oversizeCells;
    
    XML << range
	<< magnet::xml::endtag("Global");
  }

  void
  GCells::addCells(double maxdiam)
  {
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

    if (getMaxSupportedInteractionLength() < maxdiam)
      M_throw() << "The system size is too small to support the range of interactions specified (i.e. the system is smaller than the interaction diameter of one particle).";

    //Find the required size of the morton array
    magnet::math::MortonNumber<3> coords(cellCount[0], cellCount[1], cellCount[2]);
    size_t sizeReq = coords.getMortonNum();

    list.resize(sizeReq); //Empty Cells created!

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
	 << "\nVector Size <N>  " << sizeReq << std::endl;
  
    //Add the particles section
    //Required so particles find the right owning cell
    Sim->dynamics->updateAllParticles();
  
    ////Add all the particles 
    for (const size_t& id : *range)
      {
	Particle& p = Sim->particles[id];
	Sim->dynamics->updateParticle(p); 
	addToCell(id);
	if (verbose)
	  {
	    std::unordered_map<size_t, size_t>::iterator it = partCellData.find(id);
	    magnet::math::MortonNumber<3> currentCell(it->second);
	    
	    magnet::math::MortonNumber<3> estCell(getCellID(Sim->particles[ID].getPosition()));
	  
	    Vector wrapped_pos = p.getPosition();
	    for (size_t n = 0; n < NDIM; ++n)
	      {
		wrapped_pos[n] -= Sim->primaryCellSize[n] *
		  lrint(wrapped_pos[n] / Sim->primaryCellSize[n]);
	      }
	    Vector origin_pos = wrapped_pos + 0.5 * Sim->primaryCellSize - cellOffset;

	    derr << "Added particle ID=" << p.getID() << " to cell <"
		 << currentCell[0].getRealValue() 
		 << "," << currentCell[1].getRealValue()
		 << "," << currentCell[2].getRealValue()
		 << ">"
		 << "\nParticle is at this distance " << Vector(p.getPosition() - calcPosition(it->second, p)).toString() << " from the cell origin"
		 << "\nParticle position  " << p.getPosition().toString()	
		 << "\nParticle wrapped distance  " << wrapped_pos.toString()	
		 << "\nParticle relative position  " << origin_pos.toString()
		 << "\nParticle cell number  " << Vector(origin_pos[0] / cellLatticeWidth[0],
							 origin_pos[1] / cellLatticeWidth[1],
							 origin_pos[2] / cellLatticeWidth[2]
							 ).toString()
		 << std::endl;
	  }
      }

    dout << "Cell loading " << float(partCellData.size()) / NCells 
	 << std::endl;
  }

  magnet::math::MortonNumber<3>
  GCells::getCellID(Vector pos) const
  {
    Sim->BCs->applyBC(pos);

    magnet::math::MortonNumber<3> retval;

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      {
	long coord = std::floor((pos[iDim] + 0.5 * Sim->primaryCellSize[iDim] - cellOffset[iDim])
			       / cellLatticeWidth[iDim]);
	coord %= long(cellCount[iDim]);
	if (coord < 0) coord += cellCount[iDim];
	retval[iDim] = coord;
      }

    return retval;
  }

  IDRangeList
  GCells::getParticleNeighbours(const magnet::math::MortonNumber<3>& particle_cell_coords) const
  {
    if (verbose)
      {
	derr 
	  << "Getting neighbours of cell " << particle_cell_coords.toString()
	  << std::endl;
      }

    magnet::math::MortonNumber<3> zero_coords;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      zero_coords[iDim] = (particle_cell_coords[iDim].getRealValue() + cellCount[iDim] - overlink)
	% cellCount[iDim];
    
    IDRangeList retval;
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
  
  IDRangeList
  GCells::getParticleNeighbours(const Particle& part) const
  {
    return getParticleNeighbours(partCellData[part.getID()]);
  }

  IDRangeList
  GCells::getParticleNeighbours(const Vector& vec) const
  {
    return getParticleNeighbours(getCellID(vec));
  }

  double 
  GCells::getMaxSupportedInteractionLength() const
  {
    double retval(HUGE_VAL);

    for (size_t i = 0; i < NDIM; ++i)
      {
	double supported_length = cellLatticeWidth[i] * overlink
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
    Vector primaryCell = calcPosition(coords);
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

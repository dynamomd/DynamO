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

#include "packer.hpp"
#include <cmath>
#include "../datatypes/vector.hpp"
#include "../simulation/particle.hpp"
#include "../base/is_simdata.hpp"
#include "cells/include.hpp"
#include <boost/foreach.hpp>
#include <boost/random/uniform_int.hpp>
#include "../simulation/particle.hpp"
#include "../schedulers/include.hpp"
#include "../schedulers/sorters/include.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/species/include.hpp"
#include "../dynamics/units/include.hpp"
#include "../dynamics/globals/include.hpp"
#include "../dynamics/interactions/include.hpp"
#include "../dynamics/ranges/include.hpp"
#include "../dynamics/BC/include.hpp"
#include "../dynamics/liouvillean/include.hpp"
#include "../dynamics/systems/ghost.hpp"
#include "../base/is_exception.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/topology/include.hpp"
#include "../base/is_ensemble.hpp"
#include "../dynamics/locals/include.hpp"
#include "../dynamics/systems/DSMCspheres.hpp"

CIPPacker::CIPPacker(po::variables_map& vm2, DYNAMO::SimData* tmp): 
  SimBase(tmp,"SysPacker", IC_blue),
  vm(vm2)
{}

po::options_description 
CIPPacker::getOptions()
{
  po::options_description retval("System Packer General Options"), 
    hiddenopts("Packing Mode Options (description of each for each mode is "
	       "given by --packer-mode-help)");
  
  retval.add_options()
    ("packer-mode,m", po::value<size_t>(), "Chooses the system to initialise")
    ("packer-mode-help,h", 
     "Outputs the possible packer modes and their options")
    ("NCells,C", po::value<unsigned long>()->default_value(7), 
     "Number of unit cells to a dimension")
    ("xcell,x", po::value<unsigned long>(), 
     "For rectlinear co-ordinates, number of unit cells in the x direction")
    ("ycell,y", po::value<unsigned long>(), 
     "For rectlinear co-ordinates, number of unit cells in the y direction")
    ("zcell,z", po::value<unsigned long>(), 
     "For rectlinear co-ordinates, number of unit cells in the z direction")
    ("rectangular-box", "This will cause the simulation box to be deformed so "
     "that the x,y,z ecells specify the aspect ratio")
    ("density,d", po::value<Iflt>()->default_value(0.5),
     "System number density (init-mode > 1)")
    ("Thermostat,T", po::value<Iflt>(), 
     "Apply/Change the Andersen thermostat and set the Ensemble to NVT")
    //("Sentinel,S", "Installs the collision sentinal to study low densities")
    ;

  hiddenopts.add_options()
    ("b1", "boolean option one")
    ("b2", "boolean option two")
    ("i1", po::value<size_t>(), "integer option one")
    ("i2", po::value<size_t>(), "integer option two")
    ("s1", po::value<std::string>(), "string option one")
    ("s2", po::value<std::string>(), "string option two")
    ("f1", po::value<Iflt>(), "Iflt option one")
    ("f2", po::value<Iflt>(), "Iflt option two")
    ("f3", po::value<Iflt>(), "Iflt option three")
    ("f4", po::value<Iflt>(), "Iflt option four")
    ("f5", po::value<Iflt>(), "Iflt option five")
    ;

  retval.add(hiddenopts);

  return retval;
}

void
CIPPacker::initialise()
{
  if (vm.count("packer-mode-help"))
    {
      I_cout() << 
	"Modes available:\n"
	"  0: Monocomponent hard spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"  1: Monocomponent square wells\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Lambda [1.5] (well width factor)\n"
	"       --f2 : Well Depth (negative for square shoulders) [1]\n"
	"  2: Random walk of an isolated attractive polymer\n"
	"       --i1 : Chain length [20]\n"
	"       --f1 : Diameter [1.6]\n"
	"       --f2 : Well width factor [1.5]\n"
	"       --f3 : Bond inner core [0.9]\n"
	"       --f4 : Bond outer well [1.1]\n"
	"       --s1 : HP sequence to use (eg 0001010) [defaults to homopolymer if unset]\n"
	"  3: Load a config and pack it, you will need to reset the interactions etc.\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Chiral fraction (0-1) [Unloaded]\n"
	"       --s1 : File to load and use as unit cell [config.out.xml.bz2]\n"
	"  4: Monocomponent (in)elastic hard spheres in LEBC (shearing)\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Inelasticity [1.0]\n"
	"  5: Walk an isolated spiral/helix\n"
	"       --i1 : Chain length [20]\n"
	"       --i2 : Ring length (atoms in one spiral turn)[9]\n"
	"       --f1 : Diameter [1.6]\n"
	"       --f2 : Well width factor [1.5]\n"
	"       --f3 : Bond inner core (>0) [0.9]\n"
	"       --f4 : Bond outer well (>0) [1.1]\n"
	"       --f5 : Tightness of the helix, 0 is max closeness (0-1) [0.05]\n"
	"  6: Monocomponent square wells confined by two walls\n"
	"  7: Ring/Linear polymer, dropped as a straight rod\n"
	"       --i1 : Chain length (number supplied is multiplied by 2, e.g. default of 10 gives a 20mer) [10]\n"
	"       --f1 : Bond inner core (>0) [1.0]\n"
	"       --f2 : Bond outer well (>0) [1.05]\n"
	"       --f3 : Well width factor, values <= 1 use a hard sphere [1.5]\n"
	"       --b1 : If set it drops a linear chain instead of a ring\n"
	"  8: Binary Hard Spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Size Ratio (B/A), must be (0,1] [0.1]\n"
	"       --f2 : Mass Ratio (B/A) [0.001]\n"
	"       --f3 : Mol Fraction of large system (A) [0.95]\n"
	"  9: Crystal pack of lines\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
        "       --f1 : Inelasticity [1.0]\n"
	"  10: Monocomponent hard spheres using DSMC interactions\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"  11: Monocomponent hard spheres sheared using DSMC interactions\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Inelasticity [1.0]\n"
	"  12: Binary hard spheres using DSMC interactions\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Picks the g(r) to use (0:BMCSL, 1:VS, 2:HC2)\n"
	"       --f1 : Size Ratio (B/A), must be (0,1] [0.1]\n"
	"       --f2 : Mass Ratio (B/A) [0.001]\n"
	"       --f3 : Mol Fraction of large system (A) [0.95]\n"
	"  13: Crystal pack of sheared lines\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
        "       --f1 : Inelasticity [1.0]\n"
	"  14: Packing of spheres and linear rods made from stiff polymers\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Number of spheres in chain\n"
        "       --f1 : Mol fraction of spheres [0.5]\n"
        "       --f2 : Rod Length [1.0]\n"
	"  15: Monocomponent hard-parallel cubes\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --b1 : If set it enables the single occupancy model\n"
	"       --b2 : If set it bounds the system with mirrors\n"
	"  16: Stepped Potential approximating a Lennard Jones Fluid\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Sets the level of overlinking in the cell lists [1]\n"
	;
      std::cout << "\n";
      exit(1);
    }

  switch (vm["packer-mode"].as<size_t>())
    {
    case 0:
      {
	//Pack of hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));
      	
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, particleDiam, 1.0, 
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 1:
      {
	//Pack of square well molecules
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
        Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));
		
	//Set up a standard simulation
	//Just a square well system
	//old scheduler
	//Sim->ptrScheduler = new CSMultList(Sim);

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));
	
	Sim->Dynamics.setUnits(new CUSW(particleDiam,1.0, Sim));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Iflt lambda = 1.5, wellDepth = 1.0;

	if (vm.count("f1"))
	  lambda = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  wellDepth = vm["f2"].as<Iflt>();

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, particleDiam, 
						      lambda, wellDepth, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 2:
      {
	//Random walk an isolated attractive homopolymer
	size_t chainlength = 20;

	if (vm.count("i1"))
	  chainlength = vm["i1"].as<size_t>();

	Iflt sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5);

	if (vm.count("f1"))
	  sigma = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  lambda = vm["f2"].as<Iflt>();

	if (vm.count("f3"))
	  sigmin = vm["f3"].as<Iflt>();

	if (vm.count("f4"))
	  sigmax = vm["f4"].as<Iflt>();
	
	//Sit the particles 95% away of max distance from each other
	//to help with seriously overlapping wells
	Iflt diamScale = 1.0 / chainlength;
	
	CURandWalk sysPack(chainlength, (sigmin + 0.95 * (sigmax - sigmin)) 
			   * diamScale, sigma * diamScale, new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<Vector  > latticeSites(sysPack.placeObjects
					   (Vector (0,0,0)));

	//Set up the system now
	Sim->ptrScheduler = new CSDumb(Sim, new CSSBoundedPQ(Sim));
	
	Sim->Dynamics.setPBC<CNullBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction
	  (new CISquareBond(Sim, sigmin * diamScale, 
			    sigmax / sigmin,
			    new C2RChain(0, latticeSites.size()-1)
			    )
				     )->setName("Bonds");

	if (vm.count("s1"))
	  {
	    //A sequence has been supplied
	    std::vector<size_t> seq;

	    seq.resize(chainlength, 0);

	    std::string stringseq = vm["s1"].as<std::string>();

	    //Transcribe the sequence
	    for (size_t i = 0; i < chainlength; ++i)	  
	      seq[i] = boost::lexical_cast<size_t>
		(stringseq[i % stringseq.size()]);
	    
	    Sim->Dynamics.addInteraction
	      (new CISWSequence(Sim, sigma * diamScale, lambda, 1.0, 
				seq, new C2RAll()))->setName("Bulk");
	    
	    CISWSequence& interaction
	      (static_cast<CISWSequence&>
	       (*(Sim->Dynamics.getInteraction("Bulk"))));
	    
	    
	    interaction.getAlphabet().at(0).at(0) = 1.0;

	    interaction.getAlphabet().at(1).at(0) = 0.5;

	    interaction.getAlphabet().at(0).at(1) = 0.5;
	  }
	else
	  Sim->Dynamics.addInteraction(new CISquareWell(Sim, sigma * diamScale,
							lambda, 1.0, 
							1.0, 
							new C2RAll()
							))->setName("Bulk");
	
	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));

	break;
      }
    case 3:
      {
	//This packs a system using a file for the unit cell, the
	//density should just be adjusted by hand
	std::string fileName("config.out.xml.bz2");

	if (vm.count("s1"))
	  fileName = vm["s1"].as<std::string>();

	Iflt diamScale = 1.0 * vm["density"].as<Iflt>() 
	  / vm["NCells"].as<unsigned long>();

	I_cout() << "Lengthscale = " << diamScale;

	CUCell* tmpPtr;
	//Use the mirror unit cell if needed

	if (vm.count("f1"))
	  tmpPtr = new CUMirror(vm["f1"].as<Iflt>(), 
				new CUFile(Vector (diamScale,diamScale,diamScale), 
					   fileName, new CUParticle()));
	else
	  tmpPtr = new CUFile(Vector (diamScale,diamScale,diamScale), fileName, new CUParticle());
	
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(tmpPtr));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));

	Sim->Dynamics.setPBC<CSPBC>();
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, diamScale, 1.0, 
						      new C2RAll()
						      ))->setName("Bulk");
	
	
	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));
	
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 4:
      {
	//FCC simple cubic pack of hard spheres with inelasticity and shearing
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	  }

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	Iflt alpha = 1.0;

	if (vm.count("f1"))
	  alpha = vm["f1"].as<Iflt>();
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCellsShearing(Sim,"SchedulerNBList"));

	if (vm.count("rectangular-box"))
	  Sim->Dynamics.setPBC<CRLEBC>();
	else
	  Sim->Dynamics.setPBC<CSLEBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, particleDiam, alpha, 
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.setUnits(new CUShear(particleDiam, Sim));
	
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), nParticles++));

	//Insert a linear profile, zero momentum then add a vel gradient
	Sim->Dynamics.zeroMomentum(Sim->vParticleList);
	BOOST_FOREACH(CParticle& part, Sim->vParticleList)
	  part.getVelocity()[0] += part.getPosition()[1] * ShearRate;

	Sim->Ensemble.reset(new DYNAMO::CENVShear(Sim));
	break;
      }
    case 5:
      {
	//Helix/spiral layout of molecules
	size_t chainlength = 20;

	if (vm.count("i1"))
	  chainlength = vm["i1"].as<size_t>();

	size_t ringlength = 9;

	if (vm.count("i2"))
	  ringlength = vm["i2"].as<size_t>();

	Iflt sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5), 
	  tightness(0.05);

	if (vm.count("f1"))
	  sigma = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  lambda = vm["f2"].as<Iflt>();

	if (vm.count("f3"))
	  sigmin = vm["f3"].as<Iflt>();

	if (vm.count("f4"))
	  sigmax = vm["f4"].as<Iflt>();

	if (vm.count("f5"))
	  tightness = vm["f5"].as<Iflt>();
	
	//Sit the particles 95% away of max distance from each other
	//to help with seriously overlapping wells
	Iflt diamScale = 1.0 / chainlength;
	
	//Space the hard spheres 2% further apart than minimum and set
	//the bonds to 2% max length to coil this as much as possible
	CUHelix sysPack(chainlength, ringlength, 
			(sigmin + tightness * (sigmax - sigmin)) * diamScale, 
			(1.0 + tightness) * sigma * diamScale, 
			new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<Vector  > latticeSites
	  (sysPack.placeObjects(Vector (0,0,0)));

	//Set up the system now
	Sim->ptrScheduler = new CSDumb(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.setPBC<CNullBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction
	  (new CISquareBond(Sim, sigmin * diamScale, sigmax / sigmin, 
			    new C2RChain(0, latticeSites.size()-1)
			    ))->setName("Bonds");

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, sigma * diamScale, 
						      lambda, 1.0, 
						      1.0, 
						      new C2RAll()
						      ))->setName("Bulk");
	
	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 6:
      {
	//FCC simple cubic pack of hard spheres and walls
	//Pack the system, determine the number of particles
	CVector<long> cells;

	cells[0] = 3;

	cells[1] = 5;

	cells[2] = 5;

	Vector  dimensions(1,1,1);

	dimensions[0] = 0.45;
	
	boost::scoped_ptr<CUCell> sysPack(new CUFCC(cells, dimensions, 
						    new CUParticle()));
	Sim->aspectRatio[0] = 0.6;

	std::vector<Vector  > latticeSites(sysPack->placeObjects
					   (Vector (0,0,0)));
      
	Iflt particleDiam = 1.0 / 10.0;
	
	//Just a square well system
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	
	//Undo the linking of scheduler cells across the x dimension
	//D_throw() << "Needs an unlinkable scheduler";
	//static_cast<CSCells*>(Sim->ptrScheduler)->addUnlinkTask(0);

	//Cut off the x periodic boundaries
	Sim->Dynamics.setPBC<CRNoXPBC>();

	Sim->Dynamics.setUnits(new CUSW(particleDiam, 1.0, Sim));
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	Vector  norm(0,0,0), origin(0,0,0);
	norm[0] = 1.0;
	origin[0] = -0.25;

	Sim->Dynamics.addLocal(new CLWall(Sim, 1.0, norm, origin,
					  "LowWall", new CRAll(Sim)));
	/*Sim->Dynamics.addGlobal(new CGWall(Sim, 1.0, norm, origin,
	  "LowWall", new CRAll(Sim)));*/

	norm[0] = -1.0;
	origin[0] = 0.25;
	Sim->Dynamics.addLocal(new CLWall(Sim, 1.0, norm, origin, 
					   "HighWall", new CRAll(Sim)));

	/*Sim->Dynamics.addGlobal(new CGWall(Sim, 1.0, norm, origin, 
	  "HighWall", new CRAll(Sim)));*/

	Iflt lambda = 1.5;

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, particleDiam, lambda,
						      1.0, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));
	
	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 7:
      {
	//Just drops a ring polymer, you should crystalize it then
	//pack it for bulk

	size_t chainlength = 10;

	if (vm.count("i1"))
	  chainlength = vm["i1"].as<size_t>();

	Iflt sigma(1.0), sigmin(1.0), sigmax(1.05), lambda(1.5);

	if (vm.count("f1"))
	  sigmin = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  sigmax = vm["f2"].as<Iflt>();

	if (vm.count("f3"))
	  lambda = vm["f3"].as<Iflt>();
	
	//10 % more than double whats needed
	Iflt diamScale = 0.5 / (sigmax * chainlength + 2 * sigma);
	
	//CUringRod sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
	//                               * diamScale, new CUParticle());

	CUringSnake sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
			  * diamScale, new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<Vector  > latticeSites
	  (sysPack.placeObjects(Vector (0,0,0)));

	//Set up the system now
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	Sim->Dynamics.setPBC<CSPBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction
	  (new CISquareBond(Sim, sigmin * diamScale, sigmax / sigmin,
			    (vm.count("b1")) 
			    ? static_cast<C2Range*>(new C2RChain(0, latticeSites.size()-1))
			    : static_cast<C2Range*>(new C2RRing(0, latticeSites.size()-1))
			    ))->setName("Bonds");

	if (lambda >= 1.0)
	  {
	    Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	    Sim->Dynamics.addInteraction(new CISquareWell(Sim, sigma * diamScale, 
							  lambda, 1.0, 
							  1.0, 
							  new C2RAll()
							  ))->setName("Bulk");
	  }
	else
	  {
	    Sim->Dynamics.setUnits(new CUElastic(diamScale, Sim));

	    Sim->Dynamics.addInteraction(new CIHardSphere(Sim, diamScale, 1.0,
							  new C2RAll()
							  ))->setName("Bulk");
	  }
	  
	
	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "Ring"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 8:
      {
	//Pack of binary hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr
	  (new CURandomise(standardPackingHelper(new CUParticle())));

	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));
      	
	Iflt molFrac = 0.01, massFrac = 0.001, sizeRatio = 0.1;

	if (vm.count("f1"))
	  sizeRatio = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  massFrac = vm["f2"].as<Iflt>();

	if (vm.count("f3"))
	  molFrac = vm["f3"].as<Iflt>();
		  
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>()
				/ latticeSites.size(), Iflt(1.0 / 3.0));
	
	//Set up a standard simulation
	//Sim->ptrScheduler = new CSMultList(Sim);
	
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));

	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	size_t nA = static_cast<size_t>(molFrac * latticeSites.size());

	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, particleDiam, 1.0, 
			    new C2RSingle(new CRRange(0, nA - 1)))
	   )->setName("AAInt");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam, 
			    1.0, 
			    new C2RPair(new CRRange(0, nA - 1),
					new CRRange(nA, latticeSites.size()-1)))
	   )->setName("ABInt");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, sizeRatio * particleDiam, 1.0, 
			    new C2RAll()))->setName("BBInt");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(0, nA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(nA, latticeSites.size()-1),
					       massFrac, "B", 0, "BBInt")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 9:
      {
	//Pack of lines
	//Pack the system, determine the number of particles
	CURandom packroutine(vm["NCells"].as<unsigned long>(), 
			     Vector (1,1,1), Sim->uniform_sampler,
			     new CUParticle());

	packroutine.initialise();
	
	std::vector<Vector  > 
	  latticeSites(packroutine.placeObjects(Vector (0,0,0)));
      	
	Sim->Dynamics.setPBC<CSPBC>();
	  
	Iflt particleDiam = pow(vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.setLiouvillean(new CLNOrientation(Sim));

	Sim->Dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	Sim->Dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));
  
  Iflt elasticity = (vm.count("f1")) ? vm["f1"].as<Iflt>() : 1.0;

	Sim->Dynamics.addInteraction(new CILines(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSSphericalTop(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
						     particleDiam * particleDiam / 12.0,
						     "Bulk")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	static_cast<CLNOrientation&>(Sim->Dynamics.Liouvillean()).initLineOrientations(1.0);

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 10:
      {
	//Pack of DSMC hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));
      	
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>()
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	//This is to stop interactions being used for these particles
	Sim->Dynamics.addInteraction
	  (new CINull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->Dynamics.addInteraction
	  (new CIHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	//Iflt chi = 1.0 / 
	//(4.0 * tij * vm["density"].as<Iflt>() * std::sqrt(PI)); 

	Iflt packfrac = vm["density"].as<Iflt>() * PI / 6.0;

	Iflt chi = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	Iflt tij = 1.0 
	  / (4.0 * std::sqrt(PI) * vm["density"].as<Iflt>() * chi);

	//No thermostat added yet
	Sim->Dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam, 
			     2.0 * tij / latticeSites.size(), chi, 1.0, 
			     "Thermostat", new CRAll(Sim), new CRAll(Sim)));

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 11:
      {
	//Pack of DSMC hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));
      	
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt alpha = 1.0;

	if (vm.count("f1"))
	  alpha = vm["f1"].as<Iflt>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>()
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	Sim->Dynamics.setUnits(new CUShear(particleDiam, Sim));
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));
	
	Sim->Dynamics.setLiouvillean(new CLSLLOD(Sim));
	
	//This is to stop interactions being used for these particles
	Sim->Dynamics.addInteraction
	  (new CINull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->Dynamics.addInteraction
	  (new CIHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	Iflt packfrac = vm["density"].as<Iflt>() * PI / 6.0;

	Iflt chi = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	//No thermostat added yet
	Sim->Dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam, 0.001, 
			     chi, alpha, "Thermostat", 
			     new CRAll(Sim), new CRAll(Sim)));

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));
	
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 12:
      {
	//Pack of DSMC hard sphere mixture
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));
      	
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt molFrac = 0.01, massFrac = 0.001, sizeRatio = 0.1;

	if (vm.count("f1"))
	  sizeRatio = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  massFrac = vm["f2"].as<Iflt>();

	if (vm.count("f3"))
	  molFrac = vm["f3"].as<Iflt>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>()
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	//This is to stop interactions being used for these particles
	Sim->Dynamics.addInteraction
	  (new CINull(Sim, new C2RAll()))->setName("Catchall");

	size_t nA = static_cast<size_t>(molFrac * latticeSites.size());

	Iflt chiAA(1.0), chiAB(1.0), chiBB(1.0);

	size_t chimode(0);
	if (vm.count("i2"))
	  chimode = vm["i2"].as<size_t>();

	Iflt xi1 = (1.0/6.0) * PI * vm["density"].as<Iflt>()
	  * (molFrac + (1.0 - molFrac)*sizeRatio );

	Iflt xi2 = (1.0/6.0) * PI * vm["density"].as<Iflt>()
	  * (molFrac + (1.0 - molFrac)*sizeRatio*sizeRatio );

	Iflt xi3 = (1.0/6.0) * PI * vm["density"].as<Iflt>()
	  * (molFrac + (1.0 - molFrac)*sizeRatio*sizeRatio*sizeRatio);
	
	switch (chimode)
	  {
	  case 0:
	    //BMCSL
	    chiAA = (1.0/(1.0-xi3))*(1.0+3.0*xi2/(2.0*(1.0-xi3))
				     +xi2 * xi2 / (2.0*(1.0-xi3)*(1.0-xi3)));

	    chiAB = (1.0/(1.0-xi3))*(1.0+3.0*xi2/(2.0*(1.0-xi3))*sizeRatio
				     /(0.5+0.5*sizeRatio)
				     + xi2 * xi2 
				     * std::pow(sizeRatio/(0.5+0.5*sizeRatio),2)
				     / (2.0*(1.0-xi3)*(1.0-xi3)));

	    chiBB = (1.0/(1.0-xi3))*(1.0+3.0*xi2/(2.0*(1.0-xi3))*sizeRatio
				   + xi2 * xi2 *sizeRatio *sizeRatio 
				     / (2.0*(1.0-xi3)*(1.0-xi3)));
	    break;
	  case 1:
	    //VS
	    chiAA = (1.0 /(1.0 - xi3)) + (3.0 - xi3 + xi3 * xi3 * 0.5) 
	      * xi2 / (2.0 * (1.0 - xi3) * (1.0 - xi3))
	      + (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3)
	      / (6.0 * std::pow(1.0 - xi3, 3.0));

	    chiAB = (1.0 /(1.0 - xi3)) + (3.0 - xi3 + xi3 * xi3 * 0.5) 
	      * xi2 * (sizeRatio / (0.5 + 0.5*sizeRatio)) 
	      / (2.0 * (1.0 - xi3) * (1.0 - xi3))
	      + (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3)
	      * (sizeRatio / (0.5 + 0.5*sizeRatio)) 
	      * (sizeRatio / (0.5 + 0.5*sizeRatio))
	      / (6.0 * std::pow(1.0 - xi3, 3.0));

	    chiBB = (1.0 /(1.0 - xi3)) + (3.0 - xi3 + xi3 * xi3 * 0.5) 
	      * xi2 * sizeRatio / (2.0 * (1.0 - xi3) * (1.0 - xi3))
	      + (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3)
	      * sizeRatio * sizeRatio / (6.0 * std::pow(1.0 - xi3, 3.0));
	    break;
	  case 2:
	    {
	      //HC2
	      Iflt x(3.0 * (xi2 - xi3) * 0.5);

	      Iflt R(1.0 / sizeRatio);

	      chiAA = (1.0 /(1.0 - xi3)) + (3.0 - xi3 + xi3 * xi3 * 0.5) 
		* xi2 / (2.0 * (1.0 - xi3) * (1.0 - xi3))
		+ (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3)
		/ (6.0 * std::pow(1.0 - xi3, 3.0))
		+ std::exp(x) - 1.0 - x - x * x * 0.5;
	      
	      chiAB = (1.0/(1.0-xi3))*(1.0+3.0*xi2/(2.0*(1.0-xi3))*sizeRatio
				       /(0.5+0.5*sizeRatio)
				       + xi2 * xi2 
				       * std::pow(sizeRatio/(0.5+0.5*sizeRatio),2)
				       / (2.0*(1.0-xi3)*(1.0-xi3)))
		+ xi2*xi2*sizeRatio*sizeRatio*(R*R-1.0)
		/(std::pow(1.0 - xi3, 3.0)*(R+1.0)*(R+1.0))
		-xi2*xi2*xi2*sizeRatio*sizeRatio*sizeRatio*(R*R*R-1.0)
		/((1.0-xi3)*(1.0-xi3)*(1.0-xi3)*(R+1.0)*(R+1.0)*(R+1.0));
	      
	      chiBB = (1.0/(1.0-xi3))*(1.0+3.0*xi2/(2.0*(1.0-xi3))*sizeRatio
				       + xi2 * xi2 *sizeRatio *sizeRatio 
				       / (2.0*(1.0-xi3)*(1.0-xi3)));
	    }
	    break;
	  default:
	    D_throw() << "Unknown mode to set the chi's";
	  }

	chiAB *= 2.0;

	Iflt tAA = std::sqrt(PI)
	  / (chiAA * 4.0 * PI * molFrac * vm["density"].as<Iflt>());

	Iflt tAB = std::sqrt(2.0 * PI * massFrac/(1.0+massFrac))
	  / (chiAB * 4.0 * PI * (1.0 - molFrac) * vm["density"].as<Iflt>()
	     * (0.5+0.5 * sizeRatio) * (0.5+0.5 * sizeRatio));
	
	Iflt tBB = std::sqrt(PI * massFrac)
	  / (chiBB * 4.0 * PI * (1.0 - molFrac) * vm["density"].as<Iflt>()
	     * sizeRatio * sizeRatio);

	//This is to provide data on the particles
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, particleDiam, 1.0, 
			    new C2RSingle(new CRRange(0, nA - 1)))
	   )->setName("AAInt");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, sizeRatio * particleDiam, 1.0, 
			    new C2RSingle(new CRRange(nA, latticeSites.size()-1)))
	   )->setName("BBInt");

	Sim->Dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam, 
			     tAA / (2.0 * nA), chiAA, 1.0, 
			     "AADSMC", new CRRange(0, nA - 1), 
			     new CRRange(0, nA - 1)));

	Sim->Dynamics.addSystem
	  (new CSDSMCSpheres(Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam, 
			     tAB / (2.0 * nA), chiAB, 1.0, 
			     "ABDSMC", new CRRange(0, nA-1), 
			     new CRRange(nA, latticeSites.size()-1)));

	Sim->Dynamics.addSystem
	  (new CSDSMCSpheres(Sim, sizeRatio * particleDiam, 
			     tBB / (2.0 * (latticeSites.size() - nA)), chiBB, 1.0, 
			     "BBDSMC", new CRRange(nA, latticeSites.size()-1),
			     new CRRange(nA, latticeSites.size()-1)));
	

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(0, nA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(nA, latticeSites.size()-1),
					       massFrac, "B", 0, "BBInt")));

	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
      case 13:
      {
	//Pack of lines
	//Pack the system, determine the number of particles
	CURandom packroutine(vm["NCells"].as<unsigned long>(), 
			     Vector (1,1,1), Sim->uniform_sampler,
			     new CUParticle());

	packroutine.initialise();
	
	std::vector<Vector  > 
	  latticeSites(packroutine.placeObjects(Vector (0,0,0)));
      	
	Sim->Dynamics.setPBC<CSLEBC>();
	  
	Iflt particleDiam = pow(vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.setLiouvillean(new CLNOrientation(Sim));

	Sim->Dynamics.addGlobal(new CGCellsShearing(Sim,"SchedulerNBList"));
  
  Iflt elasticity = (vm.count("f1")) ? vm["f1"].as<Iflt>() : 1.0;

	Sim->Dynamics.addInteraction(new CILines(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSSphericalTop(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
						     particleDiam * particleDiam / 12.0,
						     "Bulk")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	static_cast<CLNOrientation&>(Sim->Dynamics.Liouvillean()).initLineOrientations(1.0);

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 14:
      {
	//Pack of Mings system
	//Pack the system, determine the number of particles
	Iflt molfrac(0.5), massFrac(1.0), rodlength(1.0);
	size_t chainlength(10);
	size_t nPart;

	if (vm.count("f1"))
	  molfrac = vm["f1"].as<Iflt>();

	if (vm.count("f2"))
	  rodlength = vm["f2"].as<Iflt>();

	if (vm.count("i2"))
	  chainlength = vm["i2"].as<size_t>();

	{
	  boost::scoped_ptr<CUCell> packptr
	    (standardPackingHelper(new CUParticle()));
	  
	  packptr->initialise();
	  
	  std::vector<Vector  > 
	    latticeSites(packptr->placeObjects(Vector (0,0,0)));
	  
	  nPart = latticeSites.size();
	}
	
	size_t nPartA = size_t(nPart * molfrac);

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();
	
	Iflt simVol = 1.0;
	
	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>() 
				/ nPart, 1.0 / 3.0);
	
	Iflt particleDiamB = rodlength * particleDiam / chainlength;

	boost::scoped_ptr<CUCell> packptr
	  (standardPackingHelper
	   (new CUBinary(nPartA, new CUParticle(), 
			 new CUlinearRod(chainlength, 1.05 * particleDiamB, 
					 new CUParticle()))));

	packptr->initialise();
	
	std::vector<Vector> 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));
      	
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	Sim->Dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, particleDiam, 1.0,
			    new C2RSingle(new CRRange(0, nPartA - 1)))
	   )->setName("AAInt");
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, (particleDiam + particleDiamB) / 2.0, 
			    1.0, 
			    new C2RPair
			    (new CRRange(0, nPartA - 1),
			     new CRRange(nPartA, latticeSites.size()-1)))
	   )->setName("ABInt");

	Sim->Dynamics.addInteraction
	  (new CISquareBond(Sim, 0.9 * particleDiamB, 1.1 / 0.9,
			    new C2RChains(nPartA, latticeSites.size() - 1, 
					  chainlength)
			    ))->setName("Bonds");
	
	
	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, (chainlength - 1) * particleDiamB, 1.0, 
			    new C2RChainEnds
			    (nPartA, latticeSites.size() - 1, 
			     chainlength)))->setName("RodEnds");

	Sim->Dynamics.addInteraction
	  (new CIHardSphere(Sim, particleDiamB, 1.0, 
			    new C2RAll()))->setName("BBInt");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(0, nPartA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRRange(nPartA, latticeSites.size()-1),
					       massFrac / chainlength, "B", 0, "BBInt")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
		     nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 15:
      {
	//Pack of hard spheres
	//Pack the system, determine the number of particles
	
	if (!vm.count("i1") || vm["i1"].as<size_t>() != 2)
	  D_throw() << "You should initialise cubes with simple cubic packing \"--i1 2\"";

	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));
      	
	if (latticeSites.size() % 2)
	  D_throw() << "To make sure the system has zero momentum and +-1 velocities, you must"
	    " use an even number of particles";

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>()
				/ latticeSites.size(), Iflt(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));
	Sim->Dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));


	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	if (vm.count("b1"))
	  Sim->Dynamics.addGlobal(new CGSOCells(Sim,"SOCells"));

	if (vm.count("b2"))
	  {
	    Sim->Dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(1,0,0), 
						 Vector(-0.5,0,0), 
						 "Wall1", new CRAll(Sim)));
	    Sim->Dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(0,1,0), 
						 Vector(0,-0.5,0), 
						 "Wall2", new CRAll(Sim)));
	    Sim->Dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(0,0,1), 
						 Vector(0,0,-0.5), 
						 "Wall3", new CRAll(Sim)));
	  }

	//Sim->Dynamics.addInteraction(new CIRotatedParallelCubes
	//			     (Sim, particleDiam, 1.0,
	//			      Matrix(1,0,0,0,1,0,0,0,1),
	//			      new C2RAll()))->setName("Bulk");

	Sim->Dynamics.addInteraction(new CIParallelCubes
				     (Sim, particleDiam, 1.0,
				      new C2RAll()))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, 
					       "Bulk", 0, "Bulk")));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));	
	
	size_t nParticles = 0;
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, 
		     Vector(Sim->Dynamics.units().unitVelocity(), 
			    Sim->Dynamics.units().unitVelocity(), 
			    Sim->Dynamics.units().unitVelocity()), 
		     nParticles++));
	
	{
	  boost::uniform_real<Iflt> normdist(-0.5,0.5);	  
	  boost::variate_generator<DYNAMO::baseRNG&, boost::uniform_real<Iflt> >
	    unisampler(Sim->ranGenerator, normdist);
	  
	  CVector<long> tmp = getCells();
	  
	  Vector wobblespacing;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    wobblespacing[iDim] = (Sim->aspectRatio[iDim] - particleDiam * tmp[iDim]) / tmp[iDim];
	  
	  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      part.getPosition()[iDim] += unisampler() * wobblespacing[iDim];
	}
	
	{
	  boost::variate_generator
	    <DYNAMO::baseRNG&, boost::uniform_int<unsigned int> >
	    rangen(Sim->ranGenerator, 
		   boost::uniform_int<unsigned int>
		   (0, nParticles - 1));
	  
	  size_t ID(rangen());
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    for (size_t i(0); i < nParticles / 2; ++i)
	      {
		while (Sim->vParticleList[ID].getVelocity()[iDim] < 0)
		  ID = rangen();
		
		Sim->vParticleList[ID].getVelocity()[iDim]
		  = -Sim->Dynamics.units().unitVelocity();
	      }
	}
	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 16:
      {
	//Pack of Lennard Jones stepped molecules
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector  > 
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	Iflt simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
        Iflt particleDiam = pow(simVol * vm["density"].as<Iflt>() 
				/ latticeSites.size(), Iflt(1.0 / 3.0));
		
	//Set up a standard simulation
	//Just a square well system
	//old scheduler
	//Sim->ptrScheduler = new CSMultList(Sim);

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ(Sim));

	{
	  size_t overlink = 1;
	  if (vm.count("i2"))
	    overlink = vm["i2"].as<size_t>();
	  
	  Sim->Dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList", 
					      overlink));
	}

	Sim->Dynamics.setUnits(new CUSW(particleDiam,1.0, Sim));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	typedef std::pair<Iflt,Iflt> locpair;
	std::vector<locpair> diamvec;
	diamvec.push_back(std::pair<Iflt,Iflt>(2.30,-0.06));
	diamvec.push_back(std::pair<Iflt,Iflt>(1.75,-0.22));
	diamvec.push_back(std::pair<Iflt,Iflt>(1.45,-0.55));
	diamvec.push_back(std::pair<Iflt,Iflt>(1.25,-0.98));
	diamvec.push_back(std::pair<Iflt,Iflt>(1.05,-0.47));
	diamvec.push_back(std::pair<Iflt,Iflt>(1,-0.76));
	diamvec.push_back(std::pair<Iflt,Iflt>(0.95,3.81));
	diamvec.push_back(std::pair<Iflt,Iflt>(0.90,10.95));
	diamvec.push_back(std::pair<Iflt,Iflt>(0.85,27.55));
	diamvec.push_back(std::pair<Iflt,Iflt>(0.8,66.74));
	diamvec.push_back(std::pair<Iflt,Iflt>(0.75,0));

	BOOST_FOREACH(locpair& p, diamvec)
	  {
	    p.first *= Sim->Dynamics.units().unitLength();
	    p.second *= Sim->Dynamics.units().unitEnergy();
	  }

	Sim->Dynamics.addInteraction(new CIStepped(Sim,
						   diamvec,
						   new C2RAll()
						   ))->setName("Bulk");

	Sim->Dynamics.addSpecies(smrtPlugPtr<CSpecies>
				 (new CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec() * Sim->Dynamics.units().unitVelocity(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    default:
      D_throw() << "Did not recognise the packer mode you wanted";
    }
}

void 
CIPPacker::processOptions()
{
  if (vm.count("Thermostat"))
    {      
      try {
	CSystem* thermostat = Sim->Dynamics.getSystem("Thermostat").get_ptr();

	//Only one kind of thermostat so far!
	if (dynamic_cast<const CSysGhost*>(thermostat) == NULL)
	  D_throw() << "Could not upcast thermostat to Andersens";
      
	static_cast<CSysGhost*>(thermostat)->setTemperature
	  (vm["Thermostat"].as<Iflt>() * Sim->Dynamics.units().unitEnergy());
      } catch (std::exception&)
	{
	  //No thermostat added yet
	  Sim->Dynamics.addSystem
	    (new CSysGhost(Sim, 2.0, vm["Thermostat"].as<Iflt>() 
			   * Sim->Dynamics.units().unitEnergy(), "Thermostat"));
	}

      try {
	dynamic_cast<const DYNAMO::CENVT&>(*(Sim->Ensemble));
      } catch (std::bad_cast)
	{
	  Sim->Ensemble.reset(new DYNAMO::CENVT(Sim));
	}      
    }
}

Vector  
CIPPacker::getNormalisedCellDimensions()
{
  CVector<long> cells = getCells();
  size_t maxdim = 0;
  
  //Determine the biggest dimension
  for (size_t iDim = 1; iDim < NDIM; ++iDim)
    if (cells[iDim] > cells[maxdim])
      maxdim = iDim;
  
  Vector  retval;
  
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    retval[iDim] = static_cast<Iflt>(cells[iDim])
      / static_cast<Iflt>(cells[maxdim]);
  
  return retval;
}

CUCell*
CIPPacker::standardPackingHelper(CUCell* tmpPtr)
{
  CUCell* sysPack;

  Vector  boxDimensions(1,1,1);

  if (vm.count("rectangular-box"))
    {
      boxDimensions = getNormalisedCellDimensions();
    }

  if (!vm.count("i1"))
    sysPack = new CUFCC(getCells(), boxDimensions, tmpPtr);
  else
    switch (vm["i1"].as<size_t>())
      {
      case 0:
	{
	  sysPack = new CUFCC(getCells(), boxDimensions, tmpPtr);
	  break;
	}
      case 1:
	{
	  sysPack = new CUBCC(getCells(), boxDimensions, tmpPtr);
	  break;
	}
      case 2:
	{
	  sysPack = new CUSC(getCells(), boxDimensions, tmpPtr);
	  break;
	}
      default:
	D_throw() << "Not a valid packing type (--i1)";
      }
  

  return sysPack;;
}

CVector<long> 
CIPPacker::getCells()
{
  CVector<long> cells(vm["NCells"].as<unsigned long>());

  if (vm.count("xcell"))
    cells[0] = vm["xcell"].as<unsigned long>();

  if (vm.count("ycell"))
    cells[1] = vm["ycell"].as<unsigned long>();

  if (vm.count("zcell"))
    cells[2] = vm["zcell"].as<unsigned long>();

  return cells;
}

Vector  
CIPPacker::getRandVelVec()
{
  //See http://mathworld.wolfram.com/SpherePointPicking.html
  boost::normal_distribution<Iflt> normdist(0.0, (1.0 / sqrt(NDIM)));
  
  boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution<Iflt> >
    normal_sampler(Sim->ranGenerator, normdist);
  
  Vector  tmpVec;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_sampler();
  
  return tmpVec;
}


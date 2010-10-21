/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include <magnet/exception.hpp>
#include "../base/is_simdata.hpp"
#include "../dynamics/topology/include.hpp"
#include "../base/is_ensemble.hpp"
#include "../dynamics/locals/include.hpp"
#include "../dynamics/systems/DSMCspheres.hpp"
#include "../dynamics/systems/RingDSMC.hpp"
#include "../dynamics/systems/rescale.hpp"
#include <boost/tokenizer.hpp>


CIPPacker::CIPPacker(po::variables_map& vm2, DYNAMO::SimData* tmp):
  SimBase(tmp,"SysPacker", IC_blue),
  vm(vm2)
{}

bool mySortPredictate(const Vector& v1, const Vector& v2)
{
  return v1[0] > v2[0];
}

po::options_description
CIPPacker::getOptions()
{
  po::options_description retval("System Packer General Options"),
    hiddenopts("Packing Mode Options (description of each for each mode is "
	       "given by --packer-mode-help)");

  retval.add_options()
    ("packer-mode,m", po::value<size_t>(), "Chooses the system to initialise.")
    ("packer-mode-help,h",
     "Outputs the possible packer modes and their options.")
    ("NCells,C", po::value<unsigned long>()->default_value(7),
     "Number of unit cells to a dimension.")
    ("xcell,x", po::value<unsigned long>(),
     "For rectlinear co-ordinates, number of unit cells in the x direction.")
    ("ycell,y", po::value<unsigned long>(),
     "For rectlinear co-ordinates, number of unit cells in the y direction.")
    ("zcell,z", po::value<unsigned long>(),
     "For rectlinear co-ordinates, number of unit cells in the z direction.")
    ("rectangular-box", "This will cause the simulation box to be deformed so "
     "that the x,y,z ecells specify the aspect ratio.")
    ("density,d", po::value<double>()->default_value(0.5),
     "System number density (init-mode > 1).")
    ("Thermostat,T", po::value<double>(),
     "Apply/Change the Andersen thermostat and set the Ensemble to NVT.")
    //("Sentinel,S", "Installs the collision sentinal to study low densities")
    ;

  hiddenopts.add_options()
    ("b1", "boolean option one.")
    ("b2", "boolean option two.")
    ("i1", po::value<size_t>(), "integer option one.")
    ("i2", po::value<size_t>(), "integer option two.")
    ("s1", po::value<std::string>(), "string option one.")
    ("s2", po::value<std::string>(), "string option two.")
    ("f1", po::value<double>(), "double option one.")
    ("f2", po::value<double>(), "double option two.")
    ("f3", po::value<double>(), "double option three.")
    ("f4", po::value<double>(), "double option four.")
    ("f5", po::value<double>(), "double option five.")
    ("f6", po::value<double>(), "double option six.")
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
	"       --f1 : Sets the elasticity of the hard spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Adds a temperature rescale event every x events\n"
	"       --b1 : Installs the collision sentinel for low densities\n"
	"       --b2 : Forces the use of non-morton cells in square systems\n"
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
	"  6: Monocomponent hard spheres confined by two walls, aspect ratio is set by the number of cells\n"
	"       --f1 : Elasticity of the particle and wall collisions [1]\n"	
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
	"  9: Hard needle system\n"
	"       --f1 : Inelasticity [1.0]\n"
	"       --f2 : Inertia multiplicative factor [1.0]\n"
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
	"  16: Stepped Potential\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Sets the level of overlinking in the cell lists [1]\n"
	"       --s1 : Sets the form of the stepped potential, list in r1,E1:r2,E2\n"
	"  17: Monocomponent hard spheres using Ring DSMC interactions\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Sets the fraction of T(j,k) events [1/3rd] (do not use with b1/b2)\n"
	"       --b1 : Sets chi12 to 1 [BMCSL]\n"
	"       --b2 : Sets chi13 to 1 [BMCSL]\n"
	"  18: Monocomponent sheared hard spheres using Ring DSMC interactions\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Inelasticity [0.9]\n"
	"       --b1 : Sets chi12 to 1 [BMCSL]\n"
	"       --b2 : Sets chi13 to 1 [BMCSL]\n"
	"  19: Oscillating plates bounding a system\n"
	"       --b1 : Makes the particle collisions not affect the plate\n"	
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --i2 : Upper limit on the particles inserted [All]\n"
	"       --f1 : Mass ratio [1]\n"
	"       --f2 : Length in particle radii [4.5]\n"
	"       --f3 : Hertz, if the unit of time is seconds [1]\n"
	"       --f4 : Initial displacement [130]\n"
	"       --f5 : Particle-Particle inelasticity [0.88]\n"
	"       --f6 : Particle-Wall inelasticity [0.96]\n"
	"  20: Load a set of triangles and plate it with spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --s1 : File name to load the triangles from\n"
	"       --f1 : Size scale factor of the spheres when checking for overlaps with triangles [1 = no scaling]\n"
	"  21: Pack a cylinder with spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
	"       --f1 : Length over diameter of the cylinder\n"
	"  22: Infinite system with spheres falling onto a plate with gravity\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)\n"
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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	    Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	  }
	else
	  {
	    Sim->dynamics.applyBC<BCSquarePeriodic>();

	    if (vm.count("b2"))
	      Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	    else
	      Sim->dynamics.addGlobal(new CGCellsMorton(Sim,"SchedulerNBList"));
	  }

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	if (vm.count("rectangular-box") && (vm.count("i1") && vm["i1"].as<size_t>() == 2))
	  {
	    CVector<long> cells = getCells();
	    if ((cells[0] == 1) || (cells[1] == 1) || (cells[2] == 1))
	      {
		I_cerr() << "Warning! Now assuming that you're trying to set up a 2D simulation!\n"
		  "I'm going to temporarily calculate the density by the 2D definition!";
		
		size_t dimension;
		if (cells[0] == 1)
		  dimension = 0;
		if (cells[1] == 1)
		  dimension = 1;
		if (cells[2] == 1)
		  dimension = 2;

		particleDiam = std::sqrt(simVol * vm["density"].as<double>()
					 / (Sim->aspectRatio[dimension] * latticeSites.size()));		
		
		I_cout() << "I'm changing what looks like the unused box dimension (" 
			 << dimension << ") to the optimal 2D value (3 particle diameters)";

		Sim->aspectRatio[dimension] = 3.0000001 * particleDiam;
	      }
	  }

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<MinMaxHeapPList<5> >(Sim));

	if (vm.count("b1"))
	  Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	double elasticity = 1.0;

	if (vm.count("f1"))
	  elasticity =  vm["f1"].as<double>();

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	
	if (vm.count("i2"))
	  Sim->dynamics.addSystem(new CSysRescale(Sim, vm["i2"].as<size_t>(), "RescalerEvent"));

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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

        double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	//Just a square well system
	//old scheduler
	//Sim->ptrScheduler = new CSMultList(Sim);

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));
	Sim->dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));

	Sim->dynamics.setUnits(new USquareWell(particleDiam,1.0, Sim));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	double lambda = 1.5, wellDepth = 1.0;

	if (vm.count("f1"))
	  lambda = vm["f1"].as<double>();

	if (vm.count("f2"))
	  wellDepth = vm["f2"].as<double>();

	Sim->dynamics.addInteraction(new ISquareWell(Sim, particleDiam,
						      lambda, wellDepth, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;

	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 2:
      {
	//Random walk an isolated attractive homopolymer
	size_t chainlength = 20;

	if (vm.count("i1"))
	  chainlength = vm["i1"].as<size_t>();

	double sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5);

	if (vm.count("f1"))
	  sigma = vm["f1"].as<double>();

	if (vm.count("f2"))
	  lambda = vm["f2"].as<double>();

	if (vm.count("f3"))
	  sigmin = vm["f3"].as<double>();

	if (vm.count("f4"))
	  sigmax = vm["f4"].as<double>();

	//Sit the particles 95% away of max distance from each other
	//to help with seriously overlapping wells
	double diamScale = 1.0 / chainlength;

	CURandWalk sysPack(chainlength, (sigmin + 0.95 * (sigmax - sigmin))
			   * diamScale, sigma * diamScale, new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<Vector  > latticeSites(sysPack.placeObjects
					   (Vector (0,0,0)));

	//Set up the system now
	Sim->ptrScheduler = new CSDumb(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.applyBC<BCNone>();

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction
	  (new ISquareBond(Sim, sigmin * diamScale,
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
		bool has0(false), has1(false);
	    for (size_t i = 0; i < chainlength; ++i)
	    {
	      seq[i] = boost::lexical_cast<size_t>
		(stringseq[i % stringseq.size()]);
	      if (seq[i])
			has1=true;
	      else
			has0=true;
	      if (seq[i]>1) M_throw() << "Dynamod only supports 2 types of monomers, make a sample chain and edit the configuration file by hand to use more";
	    }

	    if (has1 && has0)
	      {
		Sim->dynamics.addInteraction
		  (new ISWSequence(Sim, sigma * diamScale, lambda, 1.0,
				    seq, new C2RAll()))->setName("Bulk");
		
		ISWSequence& interaction
		  (static_cast<ISWSequence&>
		   (*(Sim->dynamics.getInteraction("Bulk"))));
		interaction.getAlphabet().at(0).at(0) = 1.0;
		
		interaction.getAlphabet().at(1).at(0) = 0.5;
		
		interaction.getAlphabet().at(0).at(1) = 0.5;
	      }
	    else if (has0 && !has1)
	      Sim->dynamics.addInteraction(new ISquareWell(Sim, sigma * diamScale,
							    lambda, 1.0,
							    1.0,
							    new C2RAll()
							    ))->setName("Bulk");
	    else if (has1 && !has0)
	      Sim->dynamics.addInteraction(new IHardSphere(Sim, sigma * diamScale, 1.0,
							    new C2RAll()
							    ))->setName("Bulk");
	  }
	else
	  Sim->dynamics.addInteraction(new ISquareWell(Sim, sigma * diamScale,
							lambda, 1.0,
							1.0,
							new C2RAll()
							))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new USquareWell(diamScale, 1.0, Sim));

	Sim->dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(), nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));

	break;
      }
    case 3:
      {
	//This packs a system using a file for the unit cell, the
	//density should just be adjusted by hand
	std::string fileName("config.out.xml.bz2");

	if (vm.count("s1"))
	  fileName = vm["s1"].as<std::string>();

	CUCell* tmpPtr;

	//Figure out how many particles are in a single unit
	tmpPtr = new CUFile(Vector (1,1,1), fileName, new CUParticle());
	tmpPtr->initialise();
	size_t NUnit = tmpPtr->placeObjects(Vector(0,0,0)).size();
	delete tmpPtr;

	//Figure out how many units there are
	tmpPtr = standardPackingHelper(new CUParticle());
	tmpPtr->initialise();
	size_t NUnitSites = tmpPtr->placeObjects(Vector(0,0,0)).size();
	delete tmpPtr;

	double diamScale = pow(vm["density"].as<double>()
			     / (NUnitSites * NUnit), double(1.0 / 3.0));

	I_cout() << "Lengthscale = " << diamScale;

	//Use the mirror unit cell if needed
	if (vm.count("f1"))
	  tmpPtr = new CUMirror(vm["f1"].as<double>(),
				new CUFile(Vector (diamScale,diamScale,diamScale),
					   fileName, new CUParticle()));
	else
	  tmpPtr = new CUFile(Vector (diamScale,diamScale,diamScale), fileName, new CUParticle());

	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(tmpPtr));
	packptr->initialise();

	std::vector<Vector  >
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));
	Sim->dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));

	Sim->dynamics.applyBC<BCSquarePeriodic>();
	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, diamScale, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");


	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new USquareWell(diamScale, 1.0, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	double alpha = 1.0;

	if (vm.count("f1"))
	  alpha = vm["f1"].as<double>();

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));
	Sim->dynamics.addGlobal(new CGCellsShearing(Sim,"SchedulerNBList"));

	if (vm.count("rectangular-box"))
	  Sim->dynamics.applyBC<BCRectangularLeesEdwards>();
	else
	  Sim->dynamics.applyBC<BCSquareLeesEdwards>();

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, alpha,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new UShear(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(), nParticles++));

	//Insert a linear profile, zero momentum then add a vel gradient
	Sim->dynamics.setCOMVelocity();
	BOOST_FOREACH(Particle& part, Sim->particleList)
	  part.getVelocity()[0] += part.getPosition()[1] * UShear::ShearRate;

	Sim->ensemble.reset(new DYNAMO::CENVShear(Sim));
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

	double sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5),
	  tightness(0.05);

	if (vm.count("f1"))
	  sigma = vm["f1"].as<double>();

	if (vm.count("f2"))
	  lambda = vm["f2"].as<double>();

	if (vm.count("f3"))
	  sigmin = vm["f3"].as<double>();

	if (vm.count("f4"))
	  sigmax = vm["f4"].as<double>();

	if (vm.count("f5"))
	  tightness = vm["f5"].as<double>();

	//Sit the particles 95% away of max distance from each other
	//to help with seriously overlapping wells
	double diamScale = 1.0 / chainlength;

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
	Sim->ptrScheduler = new CSDumb(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.applyBC<BCNone>();

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction
	  (new ISquareBond(Sim, sigmin * diamScale, sigmax / sigmin,
			    new C2RChain(0, latticeSites.size()-1)
			    ))->setName("Bonds");

	Sim->dynamics.addInteraction(new ISquareWell(Sim, sigma * diamScale,
						      lambda, 1.0,
						      1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new USquareWell(diamScale, 1.0, Sim));

	Sim->dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 6:
      {
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle(), true));
	packptr->initialise();

	std::vector<Vector  >
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));

	Sim->aspectRatio = getNormalisedCellDimensions();
	//Cut off the x periodic boundaries
	Sim->dynamics.applyBC<BCSquarePeriodicExceptX>();
	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<MinMaxHeapPList<5> >(Sim));

	if (vm.count("b1"))
	  Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));
	
	double elasticity = 1;
	if (vm.count("f1"))
	  elasticity =  vm["f1"].as<double>();

	Sim->dynamics.addLocal(new CLWall(Sim, elasticity, Vector(1,0,0), 
					  Vector(-Sim->aspectRatio[0] / 2, 0, 0),
					  "LowWall", new CRAll(Sim)));
	Sim->dynamics.addLocal(new CLWall(Sim, elasticity, Vector(-1,0,0), 
					  Vector(Sim->aspectRatio[0] / 2, 0, 0),
					  "HighWall", new CRAll(Sim)));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0, "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));

	break;
      }
    case 7:
      {
	//Just drops a ring polymer, you should crystalize it then
	//pack it for bulk

	size_t chainlength = 10;

	if (vm.count("i1"))
	  chainlength = vm["i1"].as<size_t>();

	double sigma(1.0), sigmin(1.0), sigmax(1.05), lambda(1.5);

	if (vm.count("f1"))
	  sigmin = vm["f1"].as<double>();

	if (vm.count("f2"))
	  sigmax = vm["f2"].as<double>();

	if (vm.count("f3"))
	  lambda = vm["f3"].as<double>();

	//10 % more than double whats needed
	double diamScale = 0.5 / (sigmax * chainlength + 2 * sigma);

	//CUringRod sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
	//                               * diamScale, new CUParticle());

	CUringSnake sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
			  * diamScale, new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<Vector  > latticeSites
	  (sysPack.placeObjects(Vector (0,0,0)));

	//Set up the system now
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	Sim->dynamics.applyBC<BCSquarePeriodic>();

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction
	  (new ISquareBond(Sim, sigmin * diamScale, sigmax / sigmin,
			    (vm.count("b1"))
			    ? static_cast<C2Range*>(new C2RChain(0, latticeSites.size()-1))
			    : static_cast<C2Range*>(new C2RRing(0, latticeSites.size()-1))
			    ))->setName("Bonds");

	if (lambda >= 1.0)
	  {
	    Sim->dynamics.setUnits(new USquareWell(diamScale, 1.0, Sim));

	    Sim->dynamics.addInteraction(new ISquareWell(Sim, sigma * diamScale,
							  lambda, 1.0,
							  1.0,
							  new C2RAll()
							  ))->setName("Bulk");
	  }
	else
	  {
	    Sim->dynamics.setUnits(new UHardSphere(diamScale, Sim));

	    Sim->dynamics.addInteraction(new IHardSphere(Sim, diamScale, 1.0,
							  new C2RAll()
							  ))->setName("Bulk");
	  }


	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.addStructure(new CTChain(Sim, 1, "Ring"));

	Sim->dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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

	double molFrac = 0.01, massFrac = 0.001, sizeRatio = 0.1;

	if (vm.count("f1"))
	  sizeRatio = vm["f1"].as<double>();

	if (vm.count("f2"))
	  massFrac = vm["f2"].as<double>();

	if (vm.count("f3"))
	  molFrac = vm["f3"].as<double>();

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	//Sim->ptrScheduler = new CSMultList(Sim);

	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));


	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	size_t nA = static_cast<size_t>(molFrac * latticeSites.size());

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, particleDiam, 1.0,
			    new C2RSingle(new CRRange(0, nA - 1)))
	   )->setName("AAInt");

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam,
			    1.0,
			    new C2RPair(new CRRange(0, nA - 1),
					new CRRange(nA, latticeSites.size()-1)))
	   )->setName("ABInt");

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, sizeRatio * particleDiam, 1.0,
			    new C2RAll()))->setName("BBInt");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(0, nA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(nA, latticeSites.size()-1),
					       massFrac, "B", 0, "BBInt")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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

	Sim->dynamics.applyBC<BCSquarePeriodic>();

	double particleDiam = pow(vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	//We pick a scheduler algorithm based on the density of the system
	if(vm["density"].as<double>() * 8.0 >= vm["NCells"].as<unsigned long>())
	{
	  M_throw() << "Unable to simulate systems where box volume is <= (2L)^3";
	}
	else if (vm["density"].as<double>() * 30.0 > vm["NCells"].as<unsigned long>())
	{
	  // Choose the dumb scheduler if the system volume is roughly smaller than (3L)^3
	  I_cout() << "Dumb scheduler selected due to density/particle ratio";
	  Sim->ptrScheduler = new CSDumb(Sim, new CSSBoundedPQ<>(Sim));
	}
	else
	{
	  I_cout() << "Neighbour List scheduler selected";
	  Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));
	  Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	}

	Sim->dynamics.setLiouvillean(new LNOrientation(Sim));

	Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	double elasticity = (vm.count("f1")) ? vm["f1"].as<double>() : 1.0;

	Sim->dynamics.addInteraction(new ILines(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	double inertiaMultiplicativeFactor = (vm.count("f2")) ? vm["f2"].as<double>() : 1.0;

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new SpSphericalTop(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
						     (inertiaMultiplicativeFactor * particleDiam * particleDiam) / 12.0,
						     "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	static_cast<LNOrientation&>(Sim->dynamics.getLiouvillean()).initLineOrientations(1.0);

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	//This is to stop interactions being used for these particles
	Sim->dynamics.addInteraction
	  (new INull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->dynamics.addInteraction
	  (new IHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	//double chi = 1.0 /
	//(4.0 * tij * vm["density"].as<double>() * std::sqrt(M_PI));

	double packfrac = vm["density"].as<double>() * M_PI / 6.0;

	double chi = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	double tij = 1.0
	  / (4.0 * std::sqrt(M_PI) * vm["density"].as<double>() * chi);

	//No thermostat added yet
	Sim->dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam,
			     2.0 * tij / latticeSites.size(), chi, 1.0,
			     "Thermostat", new CRAll(Sim), new CRAll(Sim)));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double alpha = 1.0;

	if (vm.count("f1"))
	  alpha = vm["f1"].as<double>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	Sim->dynamics.setUnits(new UShear(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));

	Sim->dynamics.setLiouvillean(new LSLLOD(Sim));

	//This is to stop interactions being used for these particles
	Sim->dynamics.addInteraction
	  (new INull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->dynamics.addInteraction
	  (new IHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	double packfrac = vm["density"].as<double>() * M_PI / 6.0;

	double chi = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	//No thermostat added yet
	Sim->dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam, 0.001,
			     chi, alpha, "Thermostat",
			     new CRAll(Sim), new CRAll(Sim)));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double molFrac = 0.01, massFrac = 0.001, sizeRatio = 0.1;

	if (vm.count("f1"))
	  sizeRatio = vm["f1"].as<double>();

	if (vm.count("f2"))
	  massFrac = vm["f2"].as<double>();

	if (vm.count("f3"))
	  molFrac = vm["f3"].as<double>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	//This is to stop interactions being used for these particles
	Sim->dynamics.addInteraction
	  (new INull(Sim, new C2RAll()))->setName("Catchall");

	size_t nA = static_cast<size_t>(molFrac * latticeSites.size());

	double chiAA(1.0), chiAB(1.0), chiBB(1.0);

	size_t chimode(0);
	if (vm.count("i2"))
	  chimode = vm["i2"].as<size_t>();

	double xi1 = (1.0/6.0) * M_PI * vm["density"].as<double>()
	  * (molFrac + (1.0 - molFrac)*sizeRatio );

	double xi2 = (1.0/6.0) * M_PI * vm["density"].as<double>()
	  * (molFrac + (1.0 - molFrac)*sizeRatio*sizeRatio );

	double xi3 = (1.0/6.0) * M_PI * vm["density"].as<double>()
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
	      double x(3.0 * (xi2 - xi3) * 0.5);

	      double R(1.0 / sizeRatio);

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
	    M_throw() << "Unknown mode to set the chi's";
	  }

	chiAB *= 2.0;

	double tAA = std::sqrt(M_PI)
	  / (chiAA * 4.0 * M_PI * molFrac * vm["density"].as<double>());

	double tAB = std::sqrt(2.0 * M_PI * massFrac/(1.0+massFrac))
	  / (chiAB * 4.0 * M_PI * (1.0 - molFrac) * vm["density"].as<double>()
	     * (0.5+0.5 * sizeRatio) * (0.5+0.5 * sizeRatio));

	double tBB = std::sqrt(M_PI * massFrac)
	  / (chiBB * 4.0 * M_PI * (1.0 - molFrac) * vm["density"].as<double>()
	     * sizeRatio * sizeRatio);

	//This is to provide data on the particles
	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, particleDiam, 1.0,
			    new C2RSingle(new CRRange(0, nA - 1)))
	   )->setName("AAInt");

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, sizeRatio * particleDiam, 1.0,
			    new C2RSingle(new CRRange(nA, latticeSites.size()-1)))
	   )->setName("BBInt");

	Sim->dynamics.addSystem
	  (new CSDSMCSpheres(Sim, particleDiam,
			     tAA / (2.0 * nA), chiAA, 1.0,
			     "AADSMC", new CRRange(0, nA - 1),
			     new CRRange(0, nA - 1)));

	Sim->dynamics.addSystem
	  (new CSDSMCSpheres(Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam,
			     tAB / (2.0 * nA), chiAB, 1.0,
			     "ABDSMC", new CRRange(0, nA-1),
			     new CRRange(nA, latticeSites.size()-1)));

	Sim->dynamics.addSystem
	  (new CSDSMCSpheres(Sim, sizeRatio * particleDiam,
			     tBB / (2.0 * (latticeSites.size() - nA)), chiBB, 1.0,
			     "BBDSMC", new CRRange(nA, latticeSites.size()-1),
			     new CRRange(nA, latticeSites.size()-1)));


	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(0, nA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(nA, latticeSites.size()-1),
					       massFrac, "B", 0, "BBInt")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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

	Sim->dynamics.applyBC<BCSquareLeesEdwards>();

	double particleDiam = pow(vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.setLiouvillean(new LNOrientation(Sim));

	Sim->dynamics.addGlobal(new CGCellsShearing(Sim,"SchedulerNBList"));

  double elasticity = (vm.count("f1")) ? vm["f1"].as<double>() : 1.0;

	Sim->dynamics.addInteraction(new ILines(Sim, particleDiam, elasticity,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new SpSphericalTop(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
						     particleDiam * particleDiam / 12.0,
						     "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	static_cast<LNOrientation&>(Sim->dynamics.getLiouvillean()).initLineOrientations(1.0);

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 14:
      {
	//Pack of Mings system
	//Pack the system, determine the number of particles
	double molfrac(0.5), massFrac(1.0), rodlength(1.0);
	size_t chainlength(10);
	size_t nPart;

	if (vm.count("f1"))
	  molfrac = vm["f1"].as<double>();

	if (vm.count("f2"))
	  rodlength = vm["f2"].as<double>();

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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ nPart, 1.0 / 3.0);

	double particleDiamB = rodlength * particleDiam / chainlength;

	boost::scoped_ptr<CUCell> packptr
	  (standardPackingHelper
	   (new CUBinary(nPartA, new CUParticle(),
			 new CUlinearRod(chainlength, 1.05 * particleDiamB,
					 new CUParticle()))));

	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	Sim->dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList"));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, particleDiam, 1.0,
			    new C2RSingle(new CRRange(0, nPartA - 1)))
	   )->setName("AAInt");

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, (particleDiam + particleDiamB) / 2.0,
			    1.0,
			    new C2RPair
			    (new CRRange(0, nPartA - 1),
			     new CRRange(nPartA, latticeSites.size()-1)))
	   )->setName("ABInt");

	Sim->dynamics.addInteraction
	  (new ISquareBond(Sim, 0.9 * particleDiamB, 1.1 / 0.9,
			    new C2RChains(nPartA, latticeSites.size() - 1,
					  chainlength)
			    ))->setName("Bonds");


	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, (chainlength - 1) * particleDiamB, 1.0,
			    new C2RChainEnds
			    (nPartA, latticeSites.size() - 1,
			     chainlength)))->setName("RodEnds");

	Sim->dynamics.addInteraction
	  (new IHardSphere(Sim, particleDiamB, 1.0,
			    new C2RAll()))->setName("BBInt");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(0, nPartA - 1), 1.0, "A", 0,
					       "AAInt")));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRRange(nPartA, latticeSites.size()-1),
					       massFrac / chainlength, "B", 0, "BBInt")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 15:
      {
	//Pack of hard spheres
	//Pack the system, determine the number of particles

	if (!vm.count("i1") || vm["i1"].as<size_t>() != 2)
	  M_throw() << "You should initialise cubes with simple cubic packing \"--i1 2\"";

	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));

	if (latticeSites.size() % 2)
	  M_throw() << "To make sure the system has zero momentum and +-1 velocities, you must"
	    " use an even number of particles";

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));
	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));


	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	if (vm.count("b1"))
	  Sim->dynamics.addGlobal(new CGSOCells(Sim,"SOCells"));

	if (vm.count("b2"))
	  {
	    Sim->dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(1,0,0),
						 Vector(0,0,0),
						 "Wall1", new CRAll(Sim)));
	    Sim->dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(0,1,0),
						 Vector(0,0,0),
						 "Wall2", new CRAll(Sim)));
	    Sim->dynamics.addLocal(new CLDblWall(Sim, 1.0, Vector(0,0,1),
						 Vector(0,0,0),
						 "Wall3", new CRAll(Sim)));
	  }

	//Sim->dynamics.addInteraction(new IRotatedParallelCubes
	//			     (Sim, particleDiam, 1.0,
	//			      Matrix(1,0,0,0,1,0,0,0,1),
	//			      new C2RAll()))->setName("Bulk");

	Sim->dynamics.addInteraction(new IParallelCubes
				     (Sim, particleDiam, 1.0,
				      new C2RAll()))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0,
					       "Bulk", 0, "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	size_t nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position,
		     Vector(Sim->dynamics.units().unitVelocity(),
			    Sim->dynamics.units().unitVelocity(),
			    Sim->dynamics.units().unitVelocity()),
		     nParticles++));

	{
	  boost::uniform_real<double> normdist(-0.5,0.5);
	  boost::variate_generator<DYNAMO::baseRNG&, boost::uniform_real<double> >
	    unisampler(Sim->ranGenerator, normdist);

	  CVector<long> tmp = getCells();

	  Vector wobblespacing;

	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    wobblespacing[iDim] = (Sim->aspectRatio[iDim] - particleDiam * tmp[iDim]) / tmp[iDim];

	  BOOST_FOREACH(Particle& part, Sim->particleList)
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
		while (Sim->particleList[ID].getVelocity()[iDim] < 0)
		  ID = rangen();

		Sim->particleList[ID].getVelocity()[iDim]
		  = -Sim->dynamics.units().unitVelocity();
	      }
	}
	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
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
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

        double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	//Just a square well system
	//old scheduler
	//Sim->ptrScheduler = new CSMultList(Sim);

	//New scheduler and global
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	{
	  size_t overlink = 1;
	  if (vm.count("i2"))
	    overlink = vm["i2"].as<size_t>();

	  Sim->dynamics.addGlobal(new CGCells(Sim, "SchedulerNBList",
					      overlink));
	}

	Sim->dynamics.setUnits(new USquareWell(particleDiam,1.0, Sim));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	typedef std::pair<double,double> locpair;
	std::vector<locpair> diamvec;

	if (vm.count("s1"))
	  {
	    typedef boost::tokenizer<boost::char_separator<char> >
	      tokenizer;
	    
	    tokenizer steps(vm["s1"].as<std::string>(), boost::char_separator<char>(":"));
	    BOOST_FOREACH(const std::string& step, steps)
	      {
		tokenizer stepData(step, boost::char_separator<char>(","));
		tokenizer::iterator value_iter = stepData.begin();
		locpair data;
		try {
		  data.first = boost::lexical_cast<double>(*value_iter);
		  if (++value_iter == stepData.end())
		    throw std::runtime_error("No comma");
		  data.second = boost::lexical_cast<double>(*(value_iter));
		  diamvec.push_back(data);

		  if (++value_iter != stepData.end())
		    throw std::runtime_error("Too many comma's");
		} catch (std::exception& except)
		  {
		    M_throw() << "Malformed step data, \"" << step << "\"\n" << except.what();
		  }
	      }
	  }
	else
	  {
	    diamvec.push_back(locpair(2.30, -0.06));
	    diamvec.push_back(locpair(1.75, -0.22));
	    diamvec.push_back(locpair(1.45, -0.55));
	    diamvec.push_back(locpair(1.25, -0.98));
	    diamvec.push_back(locpair(1.05, -0.47));
	    diamvec.push_back(locpair(1.00,  0.76));
	    diamvec.push_back(locpair(0.95,  3.81));
	    diamvec.push_back(locpair(0.90, 10.95));
	    diamvec.push_back(locpair(0.85, 27.55));
	    diamvec.push_back(locpair(0.80, 66.74));
	    diamvec.push_back(locpair(0.75, 1e300));
	  }
	
	I_cout() << "Building stepped potential";
	double oldr = HUGE_VAL;
	BOOST_FOREACH(locpair& p, diamvec)
	  {
	    I_cout() << "Step r=" << p.first << ", E=" << p.second;
	    if (p.first > oldr)
	      M_throw() << "Steps must be in descending order! r=" << p.first
			<< " is greater than old r=" << oldr;
	    oldr = p.first;
	    p.first *= Sim->dynamics.units().unitLength();
	    p.second *= Sim->dynamics.units().unitEnergy();
	  }

	Sim->dynamics.addInteraction(new IStepped(Sim,
						   diamvec,
						   new C2RAll()
						   ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 17:
      {
	//Pack of Ring DSMC hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	//This is to stop interactions being used for these particles
	Sim->dynamics.addInteraction
	  (new INull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->dynamics.addInteraction
	  (new IHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	double packfrac = vm["density"].as<double>() * M_PI / 6.0;

	double chi12 = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	double chi13 = chi12;

	if (vm.count("b1"))
	  chi12 = 1.0;

	if (vm.count("b2"))
	  chi13 = 1.0;

	double tij = 1.0
	  / (4.0 * std::sqrt(M_PI) * vm["density"].as<double>() * chi12);

	if (vm.count("f1"))
	  {
	    double frac = vm["f1"].as<double>();
	    chi12 = 2.0*frac*chi12;
	    chi13 = 2.0*(1.0-frac)*chi13;
	  }


	//No thermostat added yet
	Sim->dynamics.addSystem
	  (new CSRingDSMC(Sim, particleDiam,
			  2.0 * tij / latticeSites.size(), chi12, chi13, 1.0,
			  "RingDSMC", new CRAll(Sim)));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 18:
      {
	//Pack of Ring DSMC hard spheres
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector (0,0,0)));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	  }
	else
	  Sim->dynamics.applyBC<BCSquarePeriodic>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	double inelasticity = 0.9;

	if (vm.count("f1"))
	  inelasticity = vm["f1"].as<double>();

	Sim->dynamics.setUnits(new UShear(particleDiam, Sim));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSSystemOnly(Sim, new CSSCBT(Sim));

	Sim->dynamics.setLiouvillean(new LSLLOD(Sim));

	//This is to stop interactions being used for these particles
	Sim->dynamics.addInteraction
	  (new INull(Sim, new C2RAll()))->setName("Catchall");

	//This is to provide data on the particles
	Sim->dynamics.addInteraction
	  (new IHardSphere
	   (Sim, particleDiam, 1.0, new C2RAll()))->setName("Bulk");

	double packfrac = vm["density"].as<double>() * M_PI / 6.0;

	double chi12 = (1.0 - 0.5 * packfrac)
	  / std::pow(1.0 - packfrac, 3);

	double chi13 = chi12;

	if (vm.count("b1"))
	  chi12 = 1.0;

	if (vm.count("b2"))
	  chi13 = 1.0;

	double tij = 1.0
	  / (4.0 * std::sqrt(M_PI) * vm["density"].as<double>() * chi12);

	//No thermostat added yet
	Sim->dynamics.addSystem
	  (new CSRingDSMC(Sim, particleDiam,
			  2.0 * tij / latticeSites.size(), chi12, chi13, inelasticity,
			  "RingDSMC", new CRAll(Sim)));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back
	  (Particle(position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		     nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 19:
      {
	double L = 4.0;
	if (vm.count("f2"))
	  L = vm["f2"].as<double>();

	L -= 1; //This is to account for centre of mass walls

	double Delta = 13.0;
	if (vm.count("f4"))
	  Delta = vm["f4"].as<double>();

	//the  2.0 * L is to give an extra half box width on each side of the sim
	double boxL = 2.0 * L + 2.0 * Delta;
	double xy = 5.2;

	xy -= 1;//Again to account for centre of mass walls

	double Aspect =  xy / boxL;
	double MassRatio = 1.0;
	double PlateInelas = 0.96;
	if (vm.count("f6"))
	  PlateInelas = vm["f6"].as<double>();

	double ParticleInelas = 0.88;
	if (vm.count("f5"))
	  ParticleInelas = vm["f5"].as<double>();
	double boundaryInelas = PlateInelas;
	double Omega0 = M_PI * 2.0;

	if (vm.count("f1"))
	  MassRatio = vm["f1"].as<double>();

	if (vm.count("f3"))
	  Omega0 *= vm["f3"].as<double>();


	//This slight exaggeration is required to stop the cells failing with walls near the edge of the simulation
	Sim->aspectRatio = Vector(1, 1.1 * Aspect, 1.1 * Aspect);

//	Vector particleArea = Vector(0.5 * (L-2.0 * Sigma) / L ,
//				     0.9 * Aspect, 0.9 * Aspect);
//	//The minus one half spaces the particles off the wall
//	Vector particleCOM = Vector(-(0.25 * (L - 2.0 * Sigma) + Delta - 0.5)/L,
//				    0, 0);

	Vector particleArea = Vector((L + 1) / boxL, (xy + 1) / boxL, 
				     (xy + 1) / boxL);

	//The system starts at a full extention, always plus 0.1 to
	//stop instant collisions
	Vector particleCOM = Vector(Delta / boxL, 0, 0);

	CUCell* sysPack;

	if (!vm.count("i1"))
	  sysPack = new CUFCC(getCells(), particleArea, new CUParticle());
	else
	  switch (vm["i1"].as<size_t>())
	    {
	    case 0:
	      {
		sysPack = new CUFCC(getCells(), particleArea, new CUParticle());
		break;
	      }
	    case 1:
	      {
		sysPack = new CUBCC(getCells(), particleArea, new CUParticle());
		break;
	      }
	    case 2:
	      {
		sysPack = new CUSC(getCells(), particleArea, new CUParticle());
		break;
	      }
	    default:
	      M_throw() << "Not a valid packing type (--i1)";
	    }

	boost::scoped_ptr<CUCell> packptr(sysPack);
	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(particleCOM));

	Sim->dynamics.applyBC<BCNone>();
	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = 1.0 / boxL;

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	//The sentinel is needed because of the high speeds of the particles!
	Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, 
						      ParticleInelas,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addLocal(new CLWall(Sim, boundaryInelas, Vector(0,0,1), 
					  Vector(0, 0, -0.5 * Aspect),
					  "Plate2", new CRAll(Sim), false));

	Sim->dynamics.addLocal(new CLWall(Sim, boundaryInelas, Vector(0,0,-1), Vector(0, 0, +0.5 * Aspect),
					  "Plate3", new CRAll(Sim), false));

	Sim->dynamics.addLocal(new CLWall(Sim, boundaryInelas, Vector(0,+1,0), 
					  Vector(0, -0.5 * Aspect, 0),
					  "Plate4", new CRAll(Sim), false));

	Sim->dynamics.addLocal(new CLWall(Sim, boundaryInelas, Vector(0,-1,0), 
					  Vector(0, +0.5 * Aspect, 0),
					  "Plate5", new CRAll(Sim), false));

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, 
					       "Bulk", 0, "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	size_t maxPart;
	if (vm.count("i2"))
	  maxPart = vm["i2"].as<size_t>();
	else
	  maxPart = latticeSites.size();

	unsigned long nParticles = 0;
	Sim->particleList.reserve(maxPart);

	std::sort(latticeSites.begin(), latticeSites.end(), mySortPredictate);

	for (size_t i(0); i < maxPart; ++i)
	  Sim->particleList.push_back
	    (Particle(latticeSites[i], 
		       getRandVelVec() * Sim->dynamics.units().unitVelocity(),
		       nParticles++));

	bool strongPlate = false;
	if (vm.count("b1"))
	  strongPlate = true;

	Sim->dynamics.addLocal
	  (new CLOscillatingPlate(Sim, Vector(0,0,0), Vector(1,0,0), Omega0,
				  0.5 * L / boxL, PlateInelas, Delta / boxL, 
				  MassRatio * nParticles, "Plate1", 
				  new CRAll(Sim), 0.0, strongPlate));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 20:
      {
	//Pack of hard spheres then check overlaps against a set of triangles
	//Pack the system, determine the number of particles

	size_t N = boost::scoped_ptr<CUCell>(standardPackingHelper(new CUParticle()))
	  ->placeObjects(Vector(0,0,0)).size();

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->dynamics.applyBC<BCRectangularPeriodic>();
	    Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	  }
	else
	  {
	    Sim->dynamics.applyBC<BCSquarePeriodic>();
	    Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));
	  }

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>() / N, double(1.0 / 3.0));

	double overlapDiameter = particleDiam;

	if (vm.count("f1"))
	  overlapDiameter *= vm["f1"].as<double>();

	if (!vm.count("s1"))
	  M_throw() << "No triangle file name specified";

	boost::scoped_ptr<CUCell> 
	  packptr(new CUTriangleIntersect(standardPackingHelper(new CUParticle()),
					  overlapDiameter, vm["s1"].as<std::string>()));

	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	if (vm.count("b1"))
	  Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  {
	    Sim->particleList.push_back
	      (Particle(position, getRandVelVec()
			 * Sim->dynamics.units().unitVelocity(),
			 nParticles++));

	    Sim->particleList.back().getPosition()[2] -= 20 * particleDiam;
	  }

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 21:
      {
	//Pack of hard spheres in a cylinder
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();

	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));
	

	double LoverD = 1;
	if (vm.count("f1"))
	  LoverD = vm["f1"].as<double>();

	Sim->aspectRatio = Vector(1,1,1);
	
	double boxlimit;
	double cylRad = 0.5;
	if (LoverD < 1)
	  {
	    //D is unity
	    
	    //Check if the cylinder limits the sim box
	    boxlimit = LoverD;
	    if ((1.0 / std::sqrt(2.0)) < LoverD)
	      boxlimit = (1.0 / std::sqrt(2.0));

	    Sim->aspectRatio[0] = LoverD;
	  }
	else
	  {
	    //L is unity
	    Sim->aspectRatio[1] = 1.0 / LoverD;
	    Sim->aspectRatio[2] = 1.0 / LoverD;

	    boxlimit = 1.0;

	    cylRad = 0.5 / LoverD;

	    if ((1.0 / (LoverD * std::sqrt(2.0))) < 1.0)
	      boxlimit = (1.0 / (LoverD * std::sqrt(2.0)));
	  }

	//Shrink the box a little more
	boxlimit *= 0.9;

	Sim->dynamics.applyBC<BCSquarePeriodicXOnly>();
	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	double particleDiam 
	  = pow(vm["density"].as<double>() / latticeSites.size(), double(1.0 / 3.0))
	  * boxlimit;

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<>(Sim));

	if (vm.count("b1"))
	  Sim->dynamics.addGlobal(new CGPBCSentinel(Sim, "PBCSentinel"));

	Sim->dynamics.addLocal(new CLCylinder(Sim, 1.0, Vector(1,0,0), 
					      Vector(0,0,0), cylRad , "Cylinder", 
					      new CRAll(Sim), true));


	Sim->dynamics.setLiouvillean(new LNewtonian(Sim));

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					       "Bulk")));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(position * boxlimit, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
						 nParticles++));

	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 22:
      {
	//Pack of hard spheres on a plate
	//Pack the system, determine the number of particles
	boost::scoped_ptr<CUCell> packptr(standardPackingHelper(new CUParticle()));
	packptr->initialise();
	
	std::vector<Vector>
	  latticeSites(packptr->placeObjects(Vector(0,0,0)));

	Sim->aspectRatio = getNormalisedCellDimensions();
	Sim->dynamics.applyBC<BCNone>();
	Sim->dynamics.addGlobal(new CGCells(Sim,"SchedulerNBList"));

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>()
				/ latticeSites.size(), double(1.0 / 3.0));

	//Set up a standard simulation
	Sim->ptrScheduler = new CSNeighbourList(Sim, new CSSBoundedPQ<MinMaxHeapPList<5> >(Sim));

	Sim->dynamics.setUnits(new UHardSphere(particleDiam, Sim));

	Sim->dynamics.setLiouvillean(new LNewtonianGravity(Sim, -Sim->dynamics.units().unitAcceleration(), 1));

	double elasticity = 1.0;

	if (vm.count("f1"))
	  elasticity =  vm["f1"].as<double>();

	Sim->dynamics.addInteraction(new IHardSphere(Sim, particleDiam, elasticity,
						     new C2RAll()
						     ))->setName("Bulk");
	
	Sim->dynamics.addSpecies(magnet::ClonePtr<Species>
				 (new Species(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					      "Bulk")));
	
	//We actually shrink our lattice length scale by 0.999 and our
	//wall spacing by 0.9995 to prevent particles being
	//initialised touching the wall and to insert the wall just
	//inside the primary image
	Sim->dynamics.addLocal(new CLWall(Sim, 1.0, Vector(0,1,0), 
					  Vector(0, - 0.9995 * 0.5 * Sim->aspectRatio[1], 0),
					  "GroundPlate", new CRAll(Sim), false));
	
	Sim->dynamics.addGlobal(new CGParabolaSentinel(Sim,"ParabolaSentinel"));

	unsigned long nParticles = 0;
	Sim->particleList.reserve(latticeSites.size());
	BOOST_FOREACH(const Vector & position, latticeSites)
	  Sim->particleList.push_back(Particle(0.999 * position, getRandVelVec() * Sim->dynamics.units().unitVelocity(),
					       nParticles++));
	
	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
	
	break;
      }
    default:
      M_throw() << "Did not recognise the packer mode you wanted";
    }

  Sim->N = Sim->particleList.size();
}

void
CIPPacker::processOptions()
{
  if (vm.count("Thermostat"))
    {
      try {
	System* thermostat = Sim->dynamics.getSystem("Thermostat").get_ptr();

	//Only one kind of thermostat so far!
	if (dynamic_cast<const CSysGhost*>(thermostat) == NULL)
	  M_throw() << "Could not upcast thermostat to Andersens";

	static_cast<CSysGhost*>(thermostat)->setTemperature
	  (vm["Thermostat"].as<double>() * Sim->dynamics.units().unitEnergy());
      } catch (std::exception&)
	{
	  //No thermostat added yet
	  Sim->dynamics.addSystem
	    (new CSysGhost(Sim, 2.0, vm["Thermostat"].as<double>()
			   * Sim->dynamics.units().unitEnergy(), "Thermostat"));
	}

      try {
	dynamic_cast<const DYNAMO::CENVT&>(*(Sim->ensemble));
      } catch (std::bad_cast)
	{
	  Sim->ensemble.reset(new DYNAMO::CENVT(Sim));
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
    retval[iDim] = static_cast<double>(cells[iDim])
      / static_cast<double>(cells[maxdim]);

  return retval;
}

CUCell*
CIPPacker::standardPackingHelper(CUCell* tmpPtr, bool forceRectangular)
{
  CUCell* sysPack;

  Vector  boxDimensions(1,1,1);

  if (vm.count("rectangular-box") || forceRectangular)
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
	M_throw() << "Not a valid packing type (--i1)";
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
  boost::normal_distribution<double> normdist(0.0, (1.0 / sqrt(NDIM)));

  boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution<double> >
    normal_sampler(Sim->ranGenerator, normdist);

  Vector  tmpVec;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_sampler();

  return tmpVec;
}


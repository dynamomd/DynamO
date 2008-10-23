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
#include "../simulation/particle.hpp"
#include "../schedulers/include.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/species/species.hpp"
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
    ("Sentinel,S", "Installs the collision sentinal to study low densities")
    ("GCells", "Installs the cellular division global event, required by some"
     " plugins.")
    ;

  hiddenopts.add_options()
    ("i1", po::value<size_t>(), "integer option one")
    ("i2", po::value<size_t>(), "integer option two")
    ("s1", po::value<std::string>(), "string option one")
    ("s2", po::value<std::string>(), "string option two")
    ("f1", po::value<double>(), "double option one")
    ("f2", po::value<double>(), "double option two")
    ("f3", po::value<double>(), "double option three")
    ("f4", po::value<double>(), "double option four")
    ("f5", po::value<double>(), "double option five")
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
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)"
	"\n"
	"  1: Monocomponent square wells\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)"
	"\n"
	"       --f1 : Lambda [1.5] (well width factor)\n"
	"       --f2 : Well Depth (negative for square shoulders) [1]\n"
	"  2: Random walk of an isolated attractive polymer\n"
	"       --i1 : Chain length [20]\n"
	"       --f1 : Diameter [1.6]\n"
	"       --f2 : Well width factor [1.5]\n"
	"       --f3 : Bond inner core [0.9]\n"
	"       --f4 : Bond outer well [1.1]\n"
	"       --s1 : HP sequence to use (eg 0001010) "
	"[defaults to homopolymer if unset]\n"
	"  3: Load a config and pack it, you will need to reset the"
	" interactions etc.\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)"
	"\n"
	"       --f1 : Chiral fraction (0-1) [Unloaded]\n"
	"       --s1 : File to load and use as unit cell [config.out.xml.bz2]\n"
	"  4: Monocomponent (in)elastic hard spheres in LEBC (shearing)\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)"
	"\n"
	"       --f1 : Inelasticity [1.0]\n"
	"  5: Walk an isolated spiral/helix\n"
	"       --i1 : Chain length [20]\n"
	"       --i2 : Ring length (atoms in one spiral turn)[9]\n"
	"       --f1 : Diameter [1.6]\n"
	"       --f2 : Well width factor [1.5]\n"
	"       --f3 : Bond inner core (>0) [0.9]\n"
	"       --f4 : Bond outer well (>0) [1.1]\n"
	"       --f5 : Tightness of the helix, 0 is max closeness (0-1) [0.05]"
	"\n"
	"  6: Monocomponent square wells confined by two walls\n"
	"  7: Ring polymer, dropped as a straight rod\n"
	"       --i1 : Chain length (number supplied is doubled, "
	"e.g. default of 10 gives a 20mer) [10]\n"
	"       --f1 : Bond inner core (>0) [1.0]\n"
	"       --f2 : Bond outer well (>0) [1.05]\n"
	"       --f3 : Well width factor [1.5]\n"
	"  8: Binary Hard Spheres\n"
	"       --i1 : Picks the packing routine to use [0] (0:FCC,1:BCC,2:SC)"
	"\n"
	"       --f1 : Size Ratio (B/A), must be (0,1] [0.1]\n"
	"       --f2 : Mass Ratio (B/A) [0.001]\n"
	"       --f3 : Mol Fraction of large system (A) [0.95]"
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
	std::vector<CVector<> > 
	  latticeSites(standardPackingHelper(new CUParticle()));
      	
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	double particleDiam = pow(simVol * vm["density"].as<double>() 
				  / latticeSites.size(), 1.0 / 3.0);

	//Set up a standard simulation
	Sim->ptrScheduler = new CSMultList(Sim);


	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, particleDiam, 1.0, 
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					  "Bulk"));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 1:
      {
	//Pack of square well molecules
	//Pack the system, determine the number of particles
	std::vector<CVector<> > latticeSites(standardPackingHelper
					     (new CUParticle()));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
        double particleDiam = pow(simVol * vm["density"].as<double>() 
				  / latticeSites.size(), 1.0/3.0);
		
	//Set up a standard simulation
	//Just a square well system
	Sim->ptrScheduler = new CSMultList(Sim);
	
	Sim->Dynamics.setUnits(new CUSW(particleDiam,1.0, Sim));
	
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	double lambda = 1.5, wellDepth = 1.0;

	if (vm.count("f1"))
	  lambda = vm["f1"].as<double>();

	if (vm.count("f2"))
	  wellDepth = vm["f2"].as<double>();

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, particleDiam, 
						      lambda, wellDepth, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0,
					  "Bulk"));
	
	unsigned long nParticles = 0;

	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
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
	std::vector<CVector<> > latticeSites(sysPack.placeObjects
					     (CVector<>(0.0)));

	//Set up the system now
	Sim->ptrScheduler = new CSFastSingle(Sim);

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
	
	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0, 
					  "Bulk"));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec(), nParticles++));

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

	double diamScale = 1.0 * vm["density"].as<double>() 
	  / vm["NCells"].as<unsigned long>();

	I_cout() << "Lengthscale = " << diamScale;

	CUCell* tmpPtr;
	//Use the mirror unit cell if needed

	if (vm.count("f1"))
	  tmpPtr = new CUMirror(vm["f1"].as<double>(), 
				new CUFile(CVector<>(diamScale), 
					   fileName, new CUParticle()));
	else
	  tmpPtr = new CUFile(CVector<>(diamScale), fileName, new CUParticle());
	
	std::vector<CVector<> > latticeSites(standardPackingHelper(tmpPtr));
	
	//Set up the system now
	Sim->ptrScheduler = new CSMultList(Sim);
	Sim->Dynamics.setPBC<CSPBC>();
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, diamScale, 1.0, 
						      new C2RAll()
						      ))->setName("Bulk");
	
	Sim->Dynamics.addSpecies
	  (CSpecies(Sim, new CRRange(0,latticeSites.size()-1), 1.0, "Bulk", 0,
		    "Bulk"));
	
	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));
	
	unsigned long nParticles = 0;
	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 4:
      {
	//FCC simple cubic pack of hard spheres with inelasticity and shearing
	//Pack the system, determine the number of particles
	std::vector<CVector<> > latticeSites
	  (standardPackingHelper(new CUParticle()));

	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	  }

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];

	double particleDiam = pow(simVol * vm["density"].as<double>() 
				  / latticeSites.size(), 1.0 / 3.0);

	double alpha = 1.0;

	if (vm.count("f1"))
	  alpha = vm["f1"].as<double>();
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSMultListShear(Sim);

	if (vm.count("rectangular-box"))
	  Sim->Dynamics.setPBC<CRLEBC>();
	else
	  Sim->Dynamics.setPBC<CSLEBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction(new CIHardSphere(Sim, particleDiam, alpha, 
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0, 
					  "Bulk"));

	Sim->Dynamics.setUnits(new CUShear(particleDiam, Sim));
	
	unsigned long nParticles = 0;
	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec(), nParticles++));

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
	std::vector<CVector<> > latticeSites
	  (sysPack.placeObjects(CVector<>(0.0)));

	//Set up the system now
	Sim->ptrScheduler = new CSFastSingle(Sim);

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
	
	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0, 
					  "Bulk"));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "HelixPolymer"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
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

	CVector<double> dimensions(1.0);

	dimensions[0] = 0.45;
	
	boost::scoped_ptr<CUCell> sysPack(new CUFCC(cells, dimensions, 
						    new CUParticle()));
	Sim->aspectRatio[0] = 0.6;

	std::vector<CVector<> > latticeSites(sysPack->placeObjects
					     (CVector<>(0.0)));
      
	double particleDiam = 1.0 / 10.0;
	
	//Just a square well system
	Sim->ptrScheduler = new CSMultList(Sim);
	
	//Undo the linking of scheduler cells across the x dimension
	static_cast<CSCells*>(Sim->ptrScheduler)->addUnlinkTask(0);

	//Cut off the x periodic boundaries
	Sim->Dynamics.setPBC<CRNoXPBC>();

	Sim->Dynamics.setUnits(new CUSW(particleDiam, 1.0, Sim));
	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));
	
	CVector<> norm(0), origin(0);
	norm[0] = 1.0;
	origin[0] = -0.25;
	Sim->Dynamics.addGlobal(new CGWall(Sim, 1.0, norm, origin,
					   "LowWall", new CRAll(Sim)));
	norm[0] = -1.0;
	origin[0] = 0.25;
	Sim->Dynamics.addGlobal(new CGWall(Sim, 1.0, norm, origin, 
					   "HighWall", new CRAll(Sim)));

	double lambda = 1.5;

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, particleDiam, lambda,
						      1.0, 1.0,
						      new C2RAll()
						      ))->setName("Bulk");

	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0, 
					  "Bulk"));
	
	unsigned long nParticles = 0;

	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back
	  (CParticle(position, getRandVelVec(), nParticles++));

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

	double sigma(1.0), sigmin(1.0), sigmax(1.05), lambda(1.5);

	if (vm.count("f1"))
	  sigmin = vm["f1"].as<double>();

	if (vm.count("f2"))
	  sigmax = vm["f2"].as<double>();

	if (vm.count("f3"))
	  lambda = vm["f3"].as<double>();
	
	double diamScale = 1.0 / (5+chainlength);
	
	CUringRod sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
			  * diamScale, new CUParticle());

	sysPack.initialise();

	//Drop them in the middle of the sim
	std::vector<CVector<> > latticeSites
	  (sysPack.placeObjects(CVector<>(0.0)));

	//Set up the system now
	Sim->ptrScheduler = new CSFastSingle(Sim);

	Sim->Dynamics.setPBC<CNullBC>();

	Sim->Dynamics.setLiouvillean(new CLNewton(Sim));

	Sim->Dynamics.addInteraction
	  (new CISquareBond(Sim, sigmin * diamScale, sigmax / sigmin,
			    new C2RRing(0, latticeSites.size()-1)
			    ))->setName("Bonds");

	Sim->Dynamics.addInteraction(new CISquareWell(Sim, sigma * diamScale, 
						      lambda, 1.0, 
						      1.0, 
						      new C2RAll()
						      ))->setName("Bulk");
	
	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRAll(Sim), 1.0, "Bulk", 0, 
					  "Bulk"));

	Sim->Dynamics.setUnits(new CUSW(diamScale, 1.0, Sim));

	Sim->Dynamics.addStructure(new CTChain(Sim, 1, "Ring"));

	Sim->Dynamics.getTopology().back()->addMolecule(new CRAll(Sim));

	unsigned long nParticles = 0;

	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
						 nParticles++));

	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
	break;
      }
    case 8:
      {
	//Pack of hard spheres
	//Pack the system, determine the number of particles
	std::vector<CVector<> > 
	  latticeSites(standardPackingHelper(new CUParticle()));
      	
	Iflt molFrac = 0.95, massFrac = 0.001, sizeRatio = 0.1;

	if (vm.count("f1"))
	  sizeRatio = vm["f1"].as<double>();

	if (vm.count("f2"))
	  massFrac = vm["f2"].as<double>();

	if (vm.count("f3"))
	  molFrac = vm["f3"].as<double>();
		  
	if (vm.count("rectangular-box"))
	  {
	    Sim->aspectRatio = getNormalisedCellDimensions();
	    Sim->Dynamics.setPBC<CRPBC>();
	  }
	else
	  Sim->Dynamics.setPBC<CSPBC>();

	double simVol = 1.0;

	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  simVol *= Sim->aspectRatio[iDim];
	
	double particleDiam = pow(simVol * vm["density"].as<double>()
				  / latticeSites.size(), 1.0 / 3.0);
	
	//Set up a standard simulation
	Sim->ptrScheduler = new CSMultList(Sim);
	
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

	Sim->Dynamics.addSpecies(CSpecies(Sim, new CRRange(0, nA - 1), 1.0, "A",
					  0, "AAInt"));

	Sim->Dynamics.addSpecies
	  (CSpecies(Sim, new CRRange(nA, latticeSites.size()-1), 
		    massFrac, "B", 0, "BBInt"));

	Sim->Dynamics.setUnits(new CUElastic(particleDiam, Sim));
      
	unsigned long nParticles = 0;
	BOOST_FOREACH(const CVector<>& position, latticeSites)
	  Sim->vParticleList.push_back(CParticle(position, getRandVelVec(), 
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
  if (vm.count("Sentinel"))
    {      
      try {
	Sim->Dynamics.getGlobal("Cells");

	I_cout() << "Cell globals are already present.";

      } catch (std::exception&)
	{
	  Sim->Dynamics.addGlobal(new CGCells(Sim));
	}
    }

  if (vm.count("Sentinel"))
    {      
      try {
	Sim->Dynamics.getGlobal("CollisionSentinel");

	I_cout() << "Sentinel is already present.";

      } catch (std::exception&)
	{
	  Sim->Dynamics.addGlobal(new CGSentinel(Sim));
	}
    }

  if (vm.count("Thermostat"))
    {      
      try {
	CSystem* thermostat = Sim->Dynamics.getSystem("Thermostat").get_ptr();

	//Only one kind of thermostat so far!
	if (dynamic_cast<const CSysGhost*>(thermostat) == NULL)
	  D_throw() << "Could not upcast thermostat to Andersens";
      
	static_cast<const CSysGhost*>(thermostat)->setTemperature
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

CVector<> 
CIPPacker::getNormalisedCellDimensions()
{
  CVector<long> cells = getCells();
  size_t maxdim = 0;
  
  //Determine the biggest dimension
  for (size_t iDim = 1; iDim < NDIM; ++iDim)
    if (cells[iDim] > cells[maxdim])
      maxdim = iDim;
  
  CVector<> retval;
  
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    retval[iDim] = static_cast<Iflt>(cells[iDim])
      / static_cast<Iflt>(cells[maxdim]);
  
  return retval;
}

std::vector<CVector<> >
CIPPacker::standardPackingHelper(CUCell* tmpPtr)
{
  boost::scoped_ptr<CUCell> sysPack;

  CVector<> boxDimensions(1.0);

  if (vm.count("rectangular-box"))
    {
      boxDimensions = getNormalisedCellDimensions();
    }

  if (!vm.count("i1"))
    sysPack.reset(new CUFCC(getCells(), boxDimensions, tmpPtr));
  else
    switch (vm["i1"].as<size_t>())
      {
      case 0:
	{
	  sysPack.reset(new CUFCC(getCells(), boxDimensions, tmpPtr));
	  break;
	}
      case 1:
	{
	  sysPack.reset(new CUBCC(getCells(), boxDimensions, tmpPtr));
	  break;
	}
      case 2:
	{
	  sysPack.reset(new CUSC(getCells(), boxDimensions, tmpPtr));
	  break;
	}
      default:
	D_throw() << "Not a valid packing type (--i1)";
      }
  
  sysPack->initialise();
  return sysPack->placeObjects(CVector<>(0.0));
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

CVector<> 
CIPPacker::getRandVelVec()
{
  //See http://mathworld.wolfram.com/SpherePointPicking.html
  boost::normal_distribution<Iflt> normdist(0.0, (1.0 / sqrt(NDIM)));
  
  boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution<Iflt> >
    normal_sampler(Sim->ranGenerator, normdist);
  
  CVector<> tmpVec;
  for (int iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_sampler();
  
  return tmpVec;
}


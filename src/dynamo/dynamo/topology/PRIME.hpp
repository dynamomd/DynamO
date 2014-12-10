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

#pragma once
#include <dynamo/topology/topology.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <vector>

namespace dynamo {
  class TPRIME: public Topology
  {
  public:
    //All data is stored within this class to keep it scoped but available to all other classes that need it

    //Sources:
    //[1] "α-Helix formation: Discontinuous molecular dynamics on an
    //intermediate-resolution protein model", Smith and Hall (2001)
    //http://onlinelibrary.wiley.com/doi/10.1002/prot.1100/full
    //
    //[2] "Solvent effects on the conformational transition of a model
    //polyalanine peptide", Nguyen, Marchut and Hall (2004)
    //http://onlinelibrary.wiley.com/doi/10.1110/ps.04701304/full
    //
    //[3] "Spontaneous Formation of Twisted Aβ16-22 Fibrils in Large-Scale
    //Molecular-Dynamics Simulations", Cheon, Chang and Hall (2011)
    //http://www.cell.com/biophysj/fulltext/S0006-3495%2811%2901018-6
    //
    //[4] "Influence of temperature on formation of perfect tau fragment
    //fibrils using PRIME20/DMD simulations", Cheon, Chang and Hall (2012)
    //http://onlinelibrary.wiley.com/doi/10.1002/pro.2141/full
    //
    //[5] "Extending the PRIME model for protein aggregation to all 20
    //amino acids", Cheon, Chang and Hall (2010)
    //http://onlinelibrary.wiley.com/doi/10.1002/prot.22817/full
    //
    //[6] "Assembly of a tetrameric α-helical bundle: Computer simulations
    //on an intermediate-resolution protein model", Smith and Hall (2001)
    //http://onlinelibrary.wiley.com/doi/10.1002/prot.1103/abstract
    
    enum PRIME_site_type {
      NH, CH, CO, A, C, D, E, F, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
      GROUP_COUNT
    };
    
    static const std::vector<std::string> PRIME_site_names;

    //Sourced from: [3] for CH, NH, CO, K, L, V, F, A and E.
    //[4] for Q, I and Y. Other values calculated from molecular weights.
    static constexpr double _PRIME_masses[] =
      {
	//NH     CH     CO     A      C      D      E      F      H      I      K
	0.999, 0.866, 1.863, 1.000, 3.133, 3.933, 4.793, 6.061, 5.400, 3.799, 4.865,
	//L      M      N      P      Q      R      S      T      V      W      Y
	3.800, 5.000, 3.866, 2.733, 4.795, 6.666, 2.066, 3.000, 2.866, 8.666, 7.126
      };
    
    //Unbonded interaction well-depths. Sourced from [5].
    //Set to 0.000 if the interaction is hard sphere.
    static constexpr double _PRIME_well_depths[] =
      {
	/*NH-X*/0.000,0.000, 0.000, 0.000,-0.15, -0.15, -0.15,  0.000,-0.15,  0.000, 0.000, 0.000,-0.15, -0.15, 0.000,-0.15,  0.000,-0.15, -0.15,  0.000, 0.000,-0.15,
	/*CH-X*/0.000,0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	/*CO-X*/0.000,0.000, 0.000, 0.000,-0.15,  0.000, 0.000, 0.000,-0.15,  0.000,-0.15,  0.000, 0.000,-0.15, 0.000,-0.15, -0.15, -0.15, -0.15,  0.000,-0.15, -0.15,
	/*A-X*/ 0.000,0.000, 0.000,-0.084,-0.139, 0.074, 0.074,-0.148, 0.074,-0.148, 0.074,-0.148,-0.148, 0.074,0.074, 0.074, 0.074, 0.074, 0.074,-0.148,-0.148,-0.148,
	/*C-X*/-0.15, 0.000,-0.15, -0.139,-0.139,-0.116,-0.116,-0.139,-0.116,-0.139,-0.116,-0.139,-0.139,-0.116,0.015,-0.116,-0.116,-0.116,-0.116,-0.139,-0.116,-0.116,
	/*D-X*/-0.15, 0.000, 0.000, 0.074,-0.116, 0.253, 0.253, 0.015,-0.086, 0.015,-0.136, 0.015, 0.015,-0.086,0.074,-0.086,-0.136,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*E-X*/-0.15, 0.000, 0.000, 0.074,-0.116, 0.253, 0.253, 0.015,-0.086, 0.015,-0.136, 0.015, 0.015,-0.086,0.074,-0.086,-0.136,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*F-X*/ 0.000,0.000, 0.000,-0.148,-0.139, 0.015, 0.015,-0.205, 0.015,-0.203, 0.015,-0.203,-0.203, 0.015,0.015, 0.015, 0.015, 0.015, 0.015,-0.203,-0.205,-0.205,
	/*H-X*/-0.15, 0.000,-0.15,  0.074,-0.116,-0.086,-0.086, 0.015,-0.080, 0.015,-0.086, 0.015,-0.116,-0.080,0.074,-0.080,-0.086,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*I-X*/ 0.000,0.000, 0.000,-0.148,-0.139, 0.015, 0.015,-0.203, 0.015,-0.2,   0.015,-0.2,  -0.2,   0.015,0.015, 0.015, 0.015, 0.015, 0.015,-0.2,  -0.203,-0.203,
	/*K-X*/ 0.000,0.000,-0.15,  0.074,-0.116,-0.136,-0.136, 0.015,-0.086, 0.015, 0.073, 0.015,-0.116,-0.086,0.074,-0.086, 0.073,-0.086,-0.086, 0.015, 0.015,-0.086,
	/*L-X*/ 0.000,0.000, 0.000,-0.148,-0.139, 0.015, 0.015,-0.203, 0.015,-0.2,   0.015,-0.2,  -0.2,   0.015,0.015, 0.015, 0.015, 0.015, 0.015,-0.2,  -0.203,-0.203,
	/*M-X*/-0.15, 0.000, 0.000,-0.148,-0.139, 0.015, 0.015,-0.203,-0.116,-0.2,  -0.116,-0.2,  -0.2,  -0.116,0.015,-0.116,-0.116,-0.116,-0.116,-0.2,  -0.210,-0.210,
	/*N-X*/-0.15, 0.000,-0.15,  0.074,-0.116,-0.086,-0.086, 0.015,-0.080, 0.015,-0.086, 0.015,-0.116,-0.080,0.074,-0.080,-0.086,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*P-X*/ 0.000,0.000, 0.000, 0.074, 0.015, 0.074, 0.074, 0.015, 0.074, 0.015, 0.074, 0.015, 0.015, 0.074,0.074, 0.074, 0.074, 0.074, 0.074, 0.015, 0.015, 0.015,
	/*Q-X*/-0.15, 0.000,-0.15,  0.074,-0.116,-0.086,-0.086, 0.015,-0.080, 0.015,-0.086, 0.015,-0.116,-0.080,0.074,-0.080,-0.086,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*R-X*/ 0.000,0.000,-0.15,  0.074,-0.116,-0.136,-0.136, 0.015,-0.086, 0.015, 0.073, 0.015,-0.116,-0.086,0.074,-0.086, 0.073,-0.086,-0.086, 0.015, 0.015,-0.086,
	/*S-X*/-0.15, 0.000,-0.15,  0.074,-0.116,-0.086,-0.086, 0.015,-0.086, 0.015,-0.086, 0.015,-0.116,-0.086,0.074,-0.086,-0.086,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*T-X*/-0.15, 0.000,-0.15,  0.074,-0.116,-0.086,-0.086, 0.015,-0.086, 0.015,-0.086, 0.015,-0.116,-0.086,0.074,-0.086,-0.086,-0.086,-0.086, 0.015,-0.086,-0.086,
	/*V-X*/ 0.000,0.000, 0.000,-0.148,-0.139, 0.015, 0.015,-0.203, 0.015,-0.2,   0.015,-0.2,  -0.2,   0.015,0.015, 0.015, 0.015, 0.015, 0.015,-0.2,  -0.203,-0.203,
	/*W-X*/ 0.000,0.000,-0.15, -0.148,-0.116,-0.086,-0.086,-0.205,-0.086,-0.203, 0.015,-0.203,-0.210,-0.086,0.015,-0.086, 0.015,-0.086,-0.086,-0.203,-0.205,-0.201,
	/*Y-X*/-0.15, 0.000,-0.15, -0.148,-0.116,-0.086,-0.086,-0.205,-0.086,-0.203,-0.086,-0.203,-0.210,-0.086,0.015,-0.086,-0.086,-0.086,-0.086,-0.203,-0.201,-0.201
      };

    //Unbonded interaction well-diameters. Largely sourced from [5]
    //For SC-BB interactions, BB value taken from [6](!!!), SC value from [5],
    //and linear mixing rule applied.
    //Set to 0.000 if the interaction is hard sphere.
    static constexpr double _PRIME_well_diameters[] =
      {
	/*NH-X*/0.0, 0.0,0.0, 0.0,5.2,5.35,5.45,0.0,5.45,0.0,0.0, 0.0,5.45,5.25,0.0,5.4,0.0,5.3,5.35,0.0,0.0,5.6,
	/*CH-X*/0.0, 0.0,0.0, 0.0,0.0,0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,
	/*CO-X*/0.0, 0.0,0.0, 0.0,5.2,0.0, 0.0, 0.0,5.45,0.0,5.55,0.0,0.0, 5.25,0.0,5.4,5.7,5.3,5.35,0.0,5.8,5.6,
	/*A-X*/ 0.0, 0.0,0.0, 5.4,5.9,5.6, 5.9, 5.9,5.5, 5.7,6.0, 5.6,5.8, 5.6, 6.2,5.8,6.1,5.9,6.2, 6.1,5.5,5.7,
	/*C-X*/ 5.2, 0.0,5.2, 5.9,6.2,6.2, 6.1, 6.4,6.2, 6.4,6.4, 6.1,6.3, 6.2, 6.0,6.1,6.3,6.3,6.1, 6.0,6.4,6.5,
	/*D-X*/ 5.35,0.0,0.0, 5.6,6.2,6.5, 6.6, 6.7,6.6, 6.5,6.3, 6.5,6.7, 6.5, 6.3,6.3,6.5,6.1,6.2, 6.3,6.9,6.9,
	/*E-X*/ 5.45,0.0,0.0, 5.9,6.1,6.6, 6.7, 6.8,6.4, 6.6,6.4, 6.4,6.4, 6.4, 6.4,6.6,6.6,6.0,6.4, 6.5,6.9,6.8,
	/*F-X*/ 0.0, 0.0,0.0, 5.9,6.4,6.7, 6.8, 6.8,6.5, 6.6,6.9, 6.6,6.5, 6.5, 6.5,6.6,6.9,6.2,6.6, 6.5,7.0,6.8,
	/*H-X*/ 5.45,0.0,5.45,5.5,6.2,6.6, 6.4, 6.5,6.7, 6.6,6.6, 6.5,6.5, 6.5, 6.3,6.6,6.9,6.3,6.3, 6.2,7.1,6.9,
	/*I-X*/ 0.0, 0.0,0.0, 5.7,6.4,6.5, 6.6, 6.6,6.6, 6.6,6.7, 6.5,6.7, 6.6, 6.4,6.6,6.7,6.4,6.4, 6.4,6.8,6.8,
	/*K-X*/ 0.0, 0.0,5.55,6.0,6.4,6.3, 6.4, 6.9,6.6, 6.7,6.9, 6.5,6.4, 6.5, 6.7,6.7,6.8,6.1,6.5, 6.6,6.5,6.7,
	/*L-X*/ 0.0, 0.0,0.0, 5.6,6.1,6.5, 6.4, 6.6,6.5, 6.5,6.5, 6.4,6.5, 6.4, 6.3,6.3,6.8,6.3,6.2, 6.2,6.9,6.7,
	/*M-X*/ 5.45,0.0,0.0, 5.8,6.3,6.7, 6.4, 6.5,6.5, 6.7,6.4, 6.5,6.7, 6.4, 6.2,6.4,6.6,6.4,6.4, 6.4,7.0,6.6,
	/*N-X*/ 5.25,0.0,5.25,5.6,6.2,6.5, 6.4, 6.5,6.5, 6.6,6.5, 6.4,6.4, 6.3, 6.2,6.4,6.6,6.2,6.3, 6.3,6.9,6.7,
	/*P-X*/ 0.0, 0.0,0.0, 6.2,6.0,6.3, 6.4, 6.5,6.3, 6.4,6.7, 6.3,6.2, 6.2, 6.5,6.5,6.8,6.1,6.6, 6.3,6.3,6.4,
	/*Q-X*/ 5.4, 0.0,5.4, 5.8,6.1,6.3, 6.6, 6.6,6.6, 6.6,6.7, 6.3,6.4, 6.4, 6.5,6.6,6.9,6.0,6.4, 6.5,6.7,6.7,
	/*R-X*/ 0.0, 0.0,5.7, 6.1,6.3,6.5, 6.6, 6.9,6.9, 6.7,6.8, 6.8,6.6, 6.6, 6.8,6.9,7.2,6.3,6.8, 6.8,6.9,7.0,
	/*S-X*/ 5.3, 0.0,5.3, 5.9,6.3,6.1, 6.0, 6.2,6.3, 6.4,6.1, 6.3,6.4, 6.2, 6.1,6.0,6.3,6.4,6.0, 6.2,6.3,6.5,
	/*T-X*/ 5.35,0.0,5.35,6.2,6.1,6.2, 6.4, 6.6,6.3, 6.4,6.5, 6.2,6.4, 6.3, 6.6,6.4,6.8,6.0,6.5, 6.4,6.5,6.4,
	/*V-X*/ 0.0, 0.0,0.0, 6.1,6.0,6.3, 6.5, 6.5,6.2, 6.4,6.6, 6.2,6.4, 6.3, 6.3,6.5,6.8,6.2,6.4, 6.3,6.6,6.5,
	/*W-X*/ 0.0, 0.0,5.8, 5.5,6.4,6.9, 6.9, 7.0,7.1, 6.8,6.5, 6.9,7.0, 6.9, 6.3,6.7,6.9,6.3,6.5, 6.6,7.4,7.0,
	/*Y-X*/ 5.6, 0.0,5.6, 5.7,6.5,6.9, 6.8, 6.8,6.9, 6.8,6.7, 6.7,6.6, 6.7, 6.4,6.7,7.0,6.5,6.4, 6.5,7.0,7.0
      };

    //These are the bead diameters.
    //SC-SC diameters sourced from [5]. BB-BB diameters from [1].
    //SC-BB diameters parametrised to produce similar high-temperature allowed regions of the Ramachandran
    //to the alanine residue.
    static constexpr double _PRIME_diameters[] =
      {
	/*NH-X*/3.3  ,3.5  ,3.65 ,3.854,4.600,4.750,5.500,5.600,5.400,4.600,6.000,4.850,5.400,4.750,4.150,5.400,6.650,4.200,4.200,4.250,6.350,6.300,
	/*CH-X*/3.5  ,3.7  ,3.85 ,4.054,4.800,4.950,5.700,5.800,5.600,4.800,6.200,5.050,5.600,4.950,4.350,5.600,6.850,4.400,4.400,4.450,6.550,6.500,
	/*CO-X*/3.65 ,3.85 ,4.0  ,4.204,4.950,5.100,5.850,5.950,5.750,4.950,6.350,5.200,5.750,5.100,4.500,5.750,7.000,4.550,4.550,4.600,6.700,6.650,
	/*A-X*/ 3.854,4.054,4.204,2.7  ,2.8  ,2.6  ,2.9  ,2.4  ,3.1  ,2.9  ,3.3  ,2.7  ,2.9  ,2.8  ,2.9  ,3.0  ,3.0  ,2.3  ,2.6  ,2.7  ,2.7  ,2.7  ,
	/*C-X*/ 4.600,4.800,4.950,2.8  ,2.1  ,3.2  ,2.7  ,3.2  ,2.8  ,3.3  ,2.7  ,3.4  ,3.4  ,3.1  ,3.0  ,3.1  ,3.3  ,2.8  ,2.7  ,2.9  ,3.3  ,2.9  ,
	/*D-X*/ 4.750,4.950,5.100,2.6  ,3.2  ,3.4  ,2.9  ,3.1  ,2.8  ,3.4  ,3.0  ,3.0  ,3.6  ,3.2  ,3.2  ,2.8  ,3.0  ,2.8  ,3.1  ,3.0  ,3.2  ,2.8  ,
	/*E-X*/ 5.500,5.700,5.850,2.9  ,2.7  ,2.9  ,3.2  ,3.3  ,3.3  ,3.2  ,3.4  ,3.3  ,3.3  ,3.1  ,3.5  ,2.9  ,3.1  ,2.9  ,3.1  ,3.1  ,3.5  ,3.3  ,
	/*F-X*/ 5.600,5.800,5.950,2.4  ,3.2  ,3.1  ,3.3  ,3.3  ,2.9  ,3.4  ,3.5  ,3.4  ,3.2  ,2.7  ,3.1  ,3.3  ,3.3  ,2.9  ,2.8  ,3.2  ,3.4  ,3.2  ,
	/*H-X*/ 5.400,5.600,5.750,3.1  ,2.8  ,2.8  ,3.3  ,2.9  ,3.4  ,3.1  ,3.4  ,3.2  ,3.6  ,3.4  ,3.7  ,3.3  ,3.5  ,2.6  ,2.9  ,3.1  ,3.2  ,3.1  ,
	/*I-X*/ 4.600,4.800,4.950,2.9  ,3.3  ,3.4  ,3.2  ,3.4  ,3.1  ,3.3  ,2.9  ,3.4  ,3.6  ,2.8  ,3.5  ,3.1  ,3.6  ,2.6  ,3.0  ,3.3  ,3.2  ,3.0  ,
	/*K-X*/ 6.000,6.200,6.350,3.3  ,2.7  ,3.0  ,3.4  ,3.5  ,3.4  ,2.9  ,3.5  ,3.5  ,3.7  ,3.2  ,3.6  ,3.4  ,3.9  ,3.0  ,3.1  ,3.1  ,3.5  ,3.5  ,
	/*L-X*/ 4.850,5.050,5.200,2.7  ,3.4  ,3.0  ,3.3  ,3.4  ,3.2  ,3.4  ,3.5  ,3.4  ,3.6  ,3.4  ,3.5  ,3.5  ,3.4  ,3.0  ,3.2  ,3.0  ,3.4  ,3.2  ,
	/*M-X*/ 5.400,5.600,5.750,2.9  ,3.4  ,3.6  ,3.3  ,3.2  ,3.6  ,3.6  ,3.7  ,3.6  ,3.7  ,3.5  ,3.7  ,3.4  ,3.7  ,3.2  ,3.6  ,3.0  ,3.2  ,3.2  ,
	/*N-X*/ 4.750,4.950,5.100,2.8  ,3.1  ,3.2  ,3.1  ,2.7  ,3.4  ,2.8  ,3.2  ,3.4  ,3.5  ,3.3  ,3.3  ,3.5  ,2.9  ,3.0  ,3.1  ,3.1  ,2.8  ,3.3  ,
	/*P-X*/ 4.150,4.350,4.500,2.9  ,3.0  ,3.2  ,3.5  ,3.1  ,3.7  ,3.5  ,3.6  ,3.5  ,3.7  ,3.3  ,3.1  ,3.6  ,3.0  ,3.2  ,2.6  ,3.3  ,3.4  ,3.3  ,
	/*Q-X*/ 5.400,5.600,5.750,3.0  ,3.1  ,2.8  ,2.9  ,3.3  ,3.3  ,3.1  ,3.4  ,3.5  ,3.4  ,3.5  ,3.6  ,3.6  ,3.6  ,2.7  ,3.3  ,3.3  ,3.4  ,3.4  ,
	/*R-X*/ 6.650,6.850,7.000,3.0  ,3.3  ,3.0  ,3.1  ,3.3  ,3.5  ,3.6  ,3.9  ,3.4  ,3.7  ,2.9  ,3.0  ,3.6  ,3.2  ,3.0  ,3.2  ,3.1  ,3.0  ,3.1  ,
	/*S-X*/ 4.200,4.400,4.550,2.3  ,2.8  ,2.8  ,2.9  ,2.9  ,2.6  ,2.6  ,3.0  ,3.0  ,3.2  ,3.0  ,3.2  ,2.7  ,3.0  ,2.5  ,2.9  ,2.8  ,2.7  ,2.9  ,
	/*T-X*/ 4.200,4.400,4.550,2.6  ,2.7  ,3.1  ,3.1  ,2.8  ,2.9  ,3.0  ,3.1  ,3.2  ,3.6  ,3.1  ,2.6  ,3.3  ,3.2  ,2.9  ,2.9  ,2.8  ,3.3  ,3.2  ,
	/*V-X*/ 4.250,4.450,4.600,2.7  ,2.9  ,3.0  ,3.1  ,3.2  ,3.1  ,3.3  ,3.1  ,3.0  ,3.0  ,3.1  ,3.3  ,3.3  ,3.1  ,2.8  ,2.8  ,3.3  ,2.9  ,3.0  ,
	/*W-X*/ 6.350,6.550,6.700,2.7  ,3.3  ,3.2  ,3.5  ,3.4  ,3.2  ,3.2  ,3.5  ,3.4  ,3.2  ,2.8  ,3.4  ,3.4  ,3.0  ,2.7  ,3.3  ,2.9  ,3.7  ,3.2  ,
	/*Y-X*/ 6.300,6.500,6.650,2.7  ,2.9  ,2.8  ,3.3  ,3.2  ,3.1  ,3.0  ,3.5  ,3.2  ,3.2  ,3.3  ,3.3  ,3.4  ,3.1  ,2.9  ,3.2  ,3.0  ,3.2  ,3.0
      };

    //[Pseudo]Bond distances from backbone sites to the respective SC site.
    //My own parameter set: SC-CH values parametrised from PDB data, and
    //SC-{CO,NH} values from further geometric considerations.
    //
    //Each set starts with 3 dummy values so that the enum-based lookup can be used.
    static constexpr double _PRIME_SC_BB_bond_lengths[] =
      {
	//NH-SC
	//NONE  NONE  NONE  A     C     D     E     F     H     I     K
	0.00, 0.00, 0.00, 2.50, 3.20, 3.34, 4.00, 4.16, 3.93, 3.19, 4.35,
	//L     M     N     P     Q     R     S     T     V     W     Y
	3.45, 3.91, 3.34, 2.79, 3.88, 4.95, 2.83, 2.84, 2.86, 4.66, 4.59,

	//CH-SC
	//NONE  NONE  NONE  A     C     D     E     F     H     I     K
	0.00, 0.00, 0.00, 1.57, 2.37, 2.53, 3.24, 3.41, 3.16, 2.36, 3.61,
	//L     M     N     P     Q     R     S     T     V     W     Y
	2.64, 3.14, 2.53, 1.91, 3.11, 4.23, 1.96, 1.97, 1.99, 3.93, 3.86,

	//CO-SC
	//NONE  NONE  NONE  A     C     D     E     F     H     I     K
	0.00, 0.00, 0.00, 2.55, 3.25, 3.39, 4.05, 4.21, 3.98, 3.24, 4.40,
	//L     M     N     P     Q     R     S     T     V     W     Y
	3.49, 3.96, 3.39, 2.84, 3.93, 4.99, 2.88, 2.89, 2.91, 4.70, 4.64
      };

    //This is the fluctuation of the bond distance allowed.
    //Old value [1] relevant for validation purposes
    //static constexpr double _PRIME_bond_tolerance = 0.02;
    //New value [2] for more realistic Ramachandran plot
    static constexpr double _PRIME_bond_tolerance = 0.02375;

    //All subsequent values taken from [1] unless otherwise noted.

    //This is a list of the bond lengths in the backbone. These only
    //apply to beads seperated by 1 backbone bond.
    //
    //We make this a symmetric tensor to simplify the lookups. The zero
    //entries should never be used.
    static constexpr double _PRIME_BB_bond_lengths[] =
      /*        NH,    CH,    CO, */
      {/*NH*/0.000, 1.460, 1.330,
       /*CH*/1.460, 0.000, 1.510,
       /*CO*/1.330, 1.510, 0.000}
      ;

    //This is a list of the pseudobond lengths in the backbone. These
    //only apply to beads seperated by TWO backbone bonds.
    //
    //We make this a symmetric tensor to simplify the lookups. The zero
    //entries should never be used.
    static constexpr double _PRIME_pseudobond_lengths[] =
      /*       NH,   CH,   CO, */
      {/*NH*/0.00, 2.41, 2.45,
       /*CH*/2.41, 0.00, 2.45,
       /*CO*/2.45, 2.45, 0.00}
      ;

    //The next two constants are for the interactions between beads
    //separated by THREE backbone bonds.
    //
    //This is the special pseudobond length between the CH-CH groups. It
    //is the only pseudobond at this distance.
    static constexpr double _PRIME_CH_CH_pseudobond_length = 3.80;
    //
    //This is the scaling factor of the bead diameters if they are
    //closer than four bonds on the same chain.
    //Value in PRIME [1]
    //static constexpr double _PRIME_near_diameter_scale_factor = 0.75;
    //Value I've found to actually work
    static constexpr double _PRIME_near_diameter_scale_factor = 0.81;

    //Backbone-to-backbone hydrogen bonding
    //Maximum distance for the CO and NH sites
    //Old value from PRIME [1]
    //static constexpr double _PRIME_HB_well_diameter = 4.20;
    //New value from PRIME20 [4]
    static constexpr double _PRIME_HB_well_diameter = 4.50;
    //
    //Minimum distance for the "auxiliary pairs"
    static constexpr double _PRIME_HB_aux_min_distances[] =
      /*       NH,   CH,   CO, */
      {/*NH*/4.74, 5.00, 0.00,
       /*CH*/5.00, 0.00, 4.86,
       /*CO*/0.00, 4.86, 4.83}
      ;

    enum BeadLocation { NH_END, MID, CO_END };

    struct BeadData {
      BeadData(PRIME_site_type t, size_t r, BeadLocation e = MID): 
	bead_type(t), residue(r), location(e) {}

      //Define sorting for storing in maps (note location does not affect the sorting order!)
      bool operator<(const BeadData& op) const { return bead_type<op.bead_type || (!(op.bead_type<bead_type) && (residue<op.residue)); }
      bool operator==(const BeadData& op) const { return (residue == op.residue) && (bead_type == op.bead_type); }

      PRIME_site_type bead_type;
      size_t residue;
      BeadLocation location;
    };
    
    typedef boost::bimaps::bimap<boost::bimaps::unordered_set_of<size_t>, boost::bimaps::unordered_set_of<BeadData> > BeadTypeMap;

    TPRIME(const magnet::xml::Node&, dynamo::Simulation*, unsigned int ID);

    virtual ~TPRIME() {}
  
    virtual void operator<<(const magnet::xml::Node&);

    BeadData getBeadInfo(size_t ID) const { 
      const auto it = _types->left.find(ID);
#ifdef DYNAMO_DEBUG
      if (it == _types->left.end())
	M_throw() << "Particle " << ID << " has no bead data for " << getName();
#endif
      return it->second;
    }
    
    size_t getBeadID(BeadData data) const { 
      const auto it = _types->right.find(data);
#ifdef DYNAMO_DEBUG
      if (it == _types->right.end())
	M_throw() << "Particle " << ID << " has no bead data for " << getName();
#endif
      return it->second;
    }

  protected:
    std::shared_ptr<BeadTypeMap> _types;
    std::vector<std::pair<size_t, std::string> > _configData;
    virtual void outputXML(magnet::xml::XmlStream&) const;
  };

  inline size_t hash_value(const dynamo::TPRIME::BeadData& x) {
    std::size_t hash_value = 0;
    boost::hash_combine(hash_value, x.residue);
    boost::hash_combine(hash_value, x.bead_type);
    return  hash_value;
  }
}

namespace std {
  template <> struct hash<dynamo::TPRIME::BeadData>
  {
    size_t operator()(const dynamo::TPRIME::BeadData& x) const {
      return dynamo::hash_value(x);
    }
  };
}

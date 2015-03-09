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
    //For SC-BB interactions, bead diameters from _PRIME_diameters multipled by 1.5.
    //Set to 0.000 if the interaction is hard sphere.
    static constexpr double _PRIME_well_diameters[] =
      {
	/*NH-X*/0.0,0.0,0.0,0.0,6.3,6.3,6.3,0.0,6.3,0.0,0.0,0.0,6.3,6.3,0.0,6.3,0.0,6.3,6.3,0.0,0.0,6.3,
	/*CH-X*/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	/*CO-X*/0.0,0.0,0.0,0.0,6.8,0.0,0.0,0.0,6.8,0.0,6.8,0.0,0.0,6.8,0.0,6.8,6.8,6.8,6.8,0.0,6.8,6.8,
	/*A-X*/ 0.0,0.0,0.0,5.4,5.9,5.6,5.9,5.9,5.5,5.7,6.0,5.6,5.8,5.6,6.2,5.8,6.1,5.9,6.2,6.1,5.5,5.7,
	/*C-X*/ 6.3,0.0,6.8,5.9,6.2,6.2,6.1,6.4,6.2,6.4,6.4,6.1,6.3,6.2,6.0,6.1,6.3,6.3,6.1,6.0,6.4,6.5,
	/*D-X*/ 6.3,0.0,0.0,5.6,6.2,6.5,6.6,6.7,6.6,6.5,6.3,6.5,6.7,6.5,6.3,6.3,6.5,6.1,6.2,6.3,6.9,6.9,
	/*E-X*/ 6.3,0.0,0.0,5.9,6.1,6.6,6.7,6.8,6.4,6.6,6.4,6.4,6.4,6.4,6.4,6.6,6.6,6.0,6.4,6.5,6.9,6.8,
	/*F-X*/ 0.0,0.0,0.0,5.9,6.4,6.7,6.8,6.8,6.5,6.6,6.9,6.6,6.5,6.5,6.5,6.6,6.9,6.2,6.6,6.5,7.0,6.8,
	/*H-X*/ 6.3,0.0,6.8,5.5,6.2,6.6,6.4,6.5,6.7,6.6,6.6,6.5,6.5,6.5,6.3,6.6,6.9,6.3,6.3,6.2,7.1,6.9,
	/*I-X*/ 0.0,0.0,0.0,5.7,6.4,6.5,6.6,6.6,6.6,6.6,6.7,6.5,6.7,6.6,6.4,6.6,6.7,6.4,6.4,6.4,6.8,6.8,
	/*K-X*/ 0.0,0.0,6.8,6.0,6.4,6.3,6.4,6.9,6.6,6.7,6.9,6.5,6.4,6.5,6.7,6.7,6.8,6.1,6.5,6.6,6.5,6.7,
	/*L-X*/ 0.0,0.0,0.0,5.6,6.1,6.5,6.4,6.6,6.5,6.5,6.5,6.4,6.5,6.4,6.3,6.3,6.8,6.3,6.2,6.2,6.9,6.7,
	/*M-X*/ 6.3,0.0,0.0,5.8,6.3,6.7,6.4,6.5,6.5,6.7,6.4,6.5,6.7,6.4,6.2,6.4,6.6,6.4,6.4,6.4,7.0,6.6,
	/*N-X*/ 6.3,0.0,6.8,5.6,6.2,6.5,6.4,6.5,6.5,6.6,6.5,6.4,6.4,6.3,6.2,6.4,6.6,6.2,6.3,6.3,6.9,6.7,
	/*P-X*/ 0.0,0.0,0.0,6.2,6.0,6.3,6.4,6.5,6.3,6.4,6.7,6.3,6.2,6.2,6.5,6.5,6.8,6.1,6.6,6.3,6.3,6.4,
	/*Q-X*/ 6.3,0.0,6.8,5.8,6.1,6.3,6.6,6.6,6.6,6.6,6.7,6.3,6.4,6.4,6.5,6.6,6.9,6.0,6.4,6.5,6.7,6.7,
	/*R-X*/ 0.0,0.0,6.8,6.1,6.3,6.5,6.6,6.9,6.9,6.7,6.8,6.8,6.6,6.6,6.8,6.9,7.2,6.3,6.8,6.8,6.9,7.0,
	/*S-X*/ 6.3,0.0,6.8,5.9,6.3,6.1,6.0,6.2,6.3,6.4,6.1,6.3,6.4,6.2,6.1,6.0,6.3,6.4,6.0,6.2,6.3,6.5,
	/*T-X*/ 6.3,0.0,6.8,6.2,6.1,6.2,6.4,6.6,6.3,6.4,6.5,6.2,6.4,6.3,6.6,6.4,6.8,6.0,6.5,6.4,6.5,6.4,
	/*V-X*/ 0.0,0.0,0.0,6.1,6.0,6.3,6.5,6.5,6.2,6.4,6.6,6.2,6.4,6.3,6.3,6.5,6.8,6.2,6.4,6.3,6.6,6.5,
	/*W-X*/ 0.0,0.0,6.8,5.5,6.4,6.9,6.9,7.0,7.1,6.8,6.5,6.9,7.0,6.9,6.3,6.7,6.9,6.3,6.5,6.6,7.4,7.0,
	/*Y-X*/ 6.3,0.0,6.8,5.7,6.5,6.9,6.8,6.8,6.9,6.8,6.7,6.7,6.6,6.7,6.4,6.7,7.0,6.5,6.4,6.5,7.0,7.0
      };

    //These are the bead diameters.
    //SC-SC diameters sourced from [5]. BB-BB diameters from [1].
    //SC-BB diameters mix SC bead size of 5.0 linearly with the BB-BB diameters.
    //to the alanine residue.
    static constexpr double _PRIME_diameters[] =
      {
	/*NH-X*/3.3  ,3.5  ,3.65 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,4.15 ,
	/*CH-X*/3.5  ,3.7  ,3.85 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,4.35 ,
	/*CO-X*/3.65 ,3.85 ,4.0  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,4.5  ,
	/*A-X*/ 4.15 ,4.35 ,4.5  ,2.7  ,2.8  ,2.6  ,2.9  ,2.4  ,3.1  ,2.9  ,3.3  ,2.7  ,2.9  ,2.8  ,2.9  ,3.0  ,3.0  ,2.3  ,2.6  ,2.7  ,2.7  ,2.7  ,
	/*C-X*/ 4.15 ,4.35 ,4.5  ,2.8  ,2.1  ,3.2  ,2.7  ,3.2  ,2.8  ,3.3  ,2.7  ,3.4  ,3.4  ,3.1  ,3.0  ,3.1  ,3.3  ,2.8  ,2.7  ,2.9  ,3.3  ,2.9  ,
	/*D-X*/ 4.15 ,4.35 ,4.5  ,2.6  ,3.2  ,3.4  ,2.9  ,3.1  ,2.8  ,3.4  ,3.0  ,3.0  ,3.6  ,3.2  ,3.2  ,2.8  ,3.0  ,2.8  ,3.1  ,3.0  ,3.2  ,2.8  ,
	/*E-X*/ 4.15 ,4.35 ,4.5  ,2.9  ,2.7  ,2.9  ,3.2  ,3.3  ,3.3  ,3.2  ,3.4  ,3.3  ,3.3  ,3.1  ,3.5  ,2.9  ,3.1  ,2.9  ,3.1  ,3.1  ,3.5  ,3.3  ,
	/*F-X*/ 4.15 ,4.35 ,4.5  ,2.4  ,3.2  ,3.1  ,3.3  ,3.3  ,2.9  ,3.4  ,3.5  ,3.4  ,3.2  ,2.7  ,3.1  ,3.3  ,3.3  ,2.9  ,2.8  ,3.2  ,3.4  ,3.2  ,
	/*H-X*/ 4.15 ,4.35 ,4.5  ,3.1  ,2.8  ,2.8  ,3.3  ,2.9  ,3.4  ,3.1  ,3.4  ,3.2  ,3.6  ,3.4  ,3.7  ,3.3  ,3.5  ,2.6  ,2.9  ,3.1  ,3.2  ,3.1  ,
	/*I-X*/ 4.15 ,4.35 ,4.5  ,2.9  ,3.3  ,3.4  ,3.2  ,3.4  ,3.1  ,3.3  ,2.9  ,3.4  ,3.6  ,2.8  ,3.5  ,3.1  ,3.6  ,2.6  ,3.0  ,3.3  ,3.2  ,3.0  ,
	/*K-X*/ 4.15 ,4.35 ,4.5  ,3.3  ,2.7  ,3.0  ,3.4  ,3.5  ,3.4  ,2.9  ,3.5  ,3.5  ,3.7  ,3.2  ,3.6  ,3.4  ,3.9  ,3.0  ,3.1  ,3.1  ,3.5  ,3.5  ,
	/*L-X*/ 4.15 ,4.35 ,4.5  ,2.7  ,3.4  ,3.0  ,3.3  ,3.4  ,3.2  ,3.4  ,3.5  ,3.4  ,3.6  ,3.4  ,3.5  ,3.5  ,3.4  ,3.0  ,3.2  ,3.0  ,3.4  ,3.2  ,
	/*M-X*/ 4.15 ,4.35 ,4.5  ,2.9  ,3.4  ,3.6  ,3.3  ,3.2  ,3.6  ,3.6  ,3.7  ,3.6  ,3.7  ,3.5  ,3.7  ,3.4  ,3.7  ,3.2  ,3.6  ,3.0  ,3.2  ,3.2  ,
	/*N-X*/ 4.15 ,4.35 ,4.5  ,2.8  ,3.1  ,3.2  ,3.1  ,2.7  ,3.4  ,2.8  ,3.2  ,3.4  ,3.5  ,3.3  ,3.3  ,3.5  ,2.9  ,3.0  ,3.1  ,3.1  ,2.8  ,3.3  ,
	/*P-X*/ 4.15 ,4.35 ,4.5  ,2.9  ,3.0  ,3.2  ,3.5  ,3.1  ,3.7  ,3.5  ,3.6  ,3.5  ,3.7  ,3.3  ,3.1  ,3.6  ,3.0  ,3.2  ,2.6  ,3.3  ,3.4  ,3.3  ,
	/*Q-X*/ 4.15 ,4.35 ,4.5  ,3.0  ,3.1  ,2.8  ,2.9  ,3.3  ,3.3  ,3.1  ,3.4  ,3.5  ,3.4  ,3.5  ,3.6  ,3.6  ,3.6  ,2.7  ,3.3  ,3.3  ,3.4  ,3.4  ,
	/*R-X*/ 4.15 ,4.35 ,4.5  ,3.0  ,3.3  ,3.0  ,3.1  ,3.3  ,3.5  ,3.6  ,3.9  ,3.4  ,3.7  ,2.9  ,3.0  ,3.6  ,3.2  ,3.0  ,3.2  ,3.1  ,3.0  ,3.1  ,
	/*S-X*/ 4.15 ,4.35 ,4.5  ,2.3  ,2.8  ,2.8  ,2.9  ,2.9  ,2.6  ,2.6  ,3.0  ,3.0  ,3.2  ,3.0  ,3.2  ,2.7  ,3.0  ,2.5  ,2.9  ,2.8  ,2.7  ,2.9  ,
	/*T-X*/ 4.15 ,4.35 ,4.5  ,2.6  ,2.7  ,3.1  ,3.1  ,2.8  ,2.9  ,3.0  ,3.1  ,3.2  ,3.6  ,3.1  ,2.6  ,3.3  ,3.2  ,2.9  ,2.9  ,2.8  ,3.3  ,3.2  ,
	/*V-X*/ 4.15 ,4.35 ,4.5  ,2.7  ,2.9  ,3.0  ,3.1  ,3.2  ,3.1  ,3.3  ,3.1  ,3.0  ,3.0  ,3.1  ,3.3  ,3.3  ,3.1  ,2.8  ,2.8  ,3.3  ,2.9  ,3.0  ,
	/*W-X*/ 4.15 ,4.35 ,4.5  ,2.7  ,3.3  ,3.2  ,3.5  ,3.4  ,3.2  ,3.2  ,3.5  ,3.4  ,3.2  ,2.8  ,3.4  ,3.4  ,3.0  ,2.7  ,3.3  ,2.9  ,3.7  ,3.2  ,
	/*Y-X*/ 4.15 ,4.35 ,4.5  ,2.7  ,2.9  ,2.8  ,3.3  ,3.2  ,3.1  ,3.0  ,3.5  ,3.2  ,3.2  ,3.3  ,3.3  ,3.4  ,3.1  ,2.9  ,3.2  ,3.0  ,3.2  ,3.0
      };

    //[Pseudo]Bond distances from backbone sites to SC sites.
    //While the original PRIME authors use differing values for each SC site, the
    //values are not publicly available and we find this to be a suitable simplification
    static constexpr double _PRIME_SC_BB_bond_lengths[] =
      {
	//NH    CH     CO
	  2.44, 1.531, 2.49
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
    static constexpr double _PRIME_3_bonds_scale_factor = 0.75;
    //Found to be necessary for alpha-helix Hbond formation
    static constexpr double _PRIME_4_bonds_scale_factor = 0.85;

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

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

#include <dynamo/topology/PRIME.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <dynamo/property.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
  /* Out-of-class definitions for the constant data arrays*/
  const std::vector<std::string> TPRIME::PRIME_site_names {"NH", "CH", "CO", "A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

  const double TPRIME::_PRIME_diameters[22*22] = {
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

  const double TPRIME::_PRIME_masses[22] = {
    //NH     CH     CO     A      C      D      E      F      H      I      K
    0.999, 0.866, 1.863, 1.000, 3.133, 3.933, 4.793, 6.061, 5.400, 3.799, 4.865,
    //L      M      N      P      Q      R      S      T      V      W      Y
    3.800, 5.000, 3.866, 2.733, 4.795, 6.666, 2.066, 3.000, 2.866, 8.666, 7.126
  };

  const double TPRIME::_PRIME_BB_bond_lengths[9] = {
    /*        NH,    CH,    CO, */
    /*NH*/0.000, 1.460, 1.330,
    /*CH*/1.460, 0.000, 1.510,
    /*CO*/1.330, 1.510, 0.000
  };
  
  const double TPRIME::_PRIME_pseudobond_lengths[9] = {
    /*    NH,   CH,   CO, */
    /*NH*/0.00, 2.41, 2.45,
    /*CH*/2.41, 0.00, 2.45,
    /*CO*/2.45, 2.45, 0.00
  };

  const double TPRIME::_PRIME_well_diameters[22*22] = {
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

  const double TPRIME::_PRIME_well_depths[22*22] = {
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

  const double TPRIME::_PRIME_SC_BB_bond_lengths[3] = {
    //NH    CH     CO
    2.44, 1.531, 2.49
  };

  const double TPRIME::_PRIME_HB_aux_min_distances[9] = {
    /*    NH,   CH,   CO, */
    /*NH*/4.74, 5.00, 0.00,
    /*CH*/5.00, 0.00, 4.86,
    /*CO*/0.00, 4.86, 4.83
  };
  
  const double TPRIME::_PRIME_bond_tolerance = 0.02375;
  const double TPRIME::_PRIME_CH_CH_pseudobond_length = 3.80;
  const double TPRIME::_PRIME_3_bonds_scale_factor = 0.75;
  const double TPRIME::_PRIME_4_bonds_scale_factor = 0.85;
  const double TPRIME::_PRIME_HB_well_diameter = 4.50;


  /*! \brief A class which stores the type of PRIME group that the particle corresponds to.
  
    This property is added automatically when you use the PRIME topology.
  */
  class PRIMEGroupProperty: public Property
  {
    std::string _name;
    std::shared_ptr<TPRIME::BeadTypeMap> _beadTypes;

  public:
    inline PRIMEGroupProperty(const std::string name, std::shared_ptr<TPRIME::BeadTypeMap> map):
      Property(Units::Mass()), _name(name), _beadTypes(map) {}
    
    inline virtual const double getProperty(size_t ID) const
    { 
      const auto it = _beadTypes->left.find(ID);
      if (it == _beadTypes->left.end())
	M_throw() << "Do not have a PRIME bead type for particle ID " << ID;

      return TPRIME::_PRIME_masses[it->second.bead_type];
    }

    inline virtual std::string getName() const 
    { return _name; }
  
    inline virtual const double getMaxValue() const
    { return *std::max_element(TPRIME::_PRIME_masses, TPRIME::_PRIME_masses + TPRIME::GROUP_COUNT); }

    inline virtual const double getMinValue() const 
    { return *std::min_element(TPRIME::_PRIME_masses, TPRIME::_PRIME_masses + TPRIME::GROUP_COUNT); }
    //! \sa Property::rescaleUnit
    inline virtual const void rescaleUnit(const Units::Dimension dim, const double rescale)
    {
      const double factor = std::pow(rescale, _units.getUnitsPower(dim));

      if (factor && factor != 1)
	M_throw() << "Can't rescale the mass of the PRIMEGroupProperty yet!";
    }

    inline void outputParticleXMLData(magnet::xml::XmlStream& XML, const size_t pID) const
    {}
  
  
  protected:
    /*! \brief Output an XML representation of the Property to the
      passed XmlStream.
    */
    virtual void outputXML(magnet::xml::XmlStream& XML) const 
    {}
  };


  TPRIME::TPRIME(const magnet::xml::Node& XML, dynamo::Simulation* Sim, unsigned int ID):
    Topology(Sim, ID),
    _types(new BeadTypeMap)
  {
    TPRIME::operator<<(XML);
    Sim->_properties.addNamedProperty(std::shared_ptr<Property>(new PRIMEGroupProperty(_name, _types)));
  }

  void 
  TPRIME::operator<<(const magnet::xml::Node& XML) 
  {
    Topology::operator<<(XML);

    size_t residue = 0;
    for (magnet::xml::Node node = XML.findNode("Molecule"); node.valid(); ++node)
      {
	//Store the loaded data for output later (the internal
	//representation is not easy to convert back)
	_configData.push_back(std::make_pair(node.getAttribute("StartID").as<size_t>(),
					     node.getAttribute("Sequence").as<std::string>()));

	//Create an internal representation which allows fast look-ups
	const size_t startID = node.getAttribute("StartID").as<size_t>();
	size_t ID = startID;
	const std::string seq =  node.getAttribute("Sequence").as<std::string>();
	for (auto it = seq.begin(); it != seq.end(); ++it)
	  {
	    BeadLocation location = MID;
	    if (it == seq.begin()) location = NH_END;
	    if (it == (seq.end() - 1)) location = CO_END;

	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(NH, residue, location)));
	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(CH, residue, location)));
	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(CO, residue, location)));
	    
	    switch (*it) {
	    case 'A': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::A, residue, location))); break;
	    case 'C': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::C, residue, location))); break;
	    case 'D': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::D, residue, location))); break;
	    case 'E': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::E, residue, location))); break;
	    case 'F': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::F, residue, location))); break;
	    case 'G': /*This residue has no side chain*/ break;
	    case 'H': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::H, residue, location))); break;
	    case 'I': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::I, residue, location))); break;
	    case 'K': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::K, residue, location))); break;
	    case 'L': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::L, residue, location))); break;
	    case 'M': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::M, residue, location))); break;
	    case 'N': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::N, residue, location))); break;
	    case 'P': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::P, residue, location))); break;
	    case 'Q': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::Q, residue, location))); break;
	    case 'R': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::R, residue, location))); break;
	    case 'S': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::S, residue, location))); break;
	    case 'T': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::T, residue, location))); break;
	    case 'V': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::V, residue, location))); break;
	    case 'W': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::W, residue, location))); break;
	    case 'Y': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::Y, residue, location))); break;
	    default:
	      M_throw() << "Unrecognised PRIME group type " << *it;
	    }
	    ++residue;
	  }

	ranges.push_back(shared_ptr<IDRange>(new IDRangeRange(startID, ID)));

	//If we're starting a new chain, skip at least three residue
	//ID's (so that all special cases for intra-molecule
	//interactions cannot happen between molecules)
	residue += 3;
      }
  }

  void 
  TPRIME::outputXML(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::attr("Name") << _name;
    XML << magnet::xml::attr("Type") << "PRIME";

    for (const auto& entry: _configData)
      XML << magnet::xml::tag("Molecule")
	  << magnet::xml::attr("StartID") << entry.first
	  << magnet::xml::attr("Sequence") << entry.second
	  << magnet::xml::endtag("Molecule");
  }
}

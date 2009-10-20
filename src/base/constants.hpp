/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <cstddef>

#ifndef Consts_h
# define Consts_h
# ifdef DYNAMO_double_precsision
   typedef double Iflt;
   typedef long double lIflt;
   const Iflt eps(1e-8);
# else
   typedef float Iflt;
   typedef float lIflt;
   const Iflt eps(1e-2);
# endif
 const Iflt PI(3.14159265358979323846);
 const size_t NDIM(3);
 const Iflt ShearRate(1);
 const Iflt Gravity(10);
 const char configFileVersion[] = "1.2.0";
#endif

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

#ifndef Consts_h
# define Consts_h
# define NDIM 3
# define PI 3.14159265358979323846
# define ShearRate 1
# define configFileVersion "1.1.3"
# ifdef DYNAMO_double_precsision
#  define eps 1e-8
   typedef double Iflt;
# else
#  define eps 1e-2
   typedef float Iflt;
# endif
#endif

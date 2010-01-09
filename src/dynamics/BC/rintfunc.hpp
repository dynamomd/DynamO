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
#pragma once

#define	FLOAT_TO_INT(in,out)		\
	__asm__ __volatile__ ("fistpl %0" : "=m" (out) : "t" (in) : "st") ;

inline int rintfunc(const Iflt& input)
{
  int retval;
  FLOAT_TO_INT(input, retval);
  return retval;
}

//#ifdef DYNAMO_double_precsision
//# define rintfunc lrint
//#else
//# define rintfunc lrintf
//#endif

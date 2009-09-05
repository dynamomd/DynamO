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

#pragma once

#include "../../../extcode/mathtemplates.hpp"

template<class T>
Iflt quadRootHunter(const T& fL, Iflt length, Iflt& t_low, Iflt& t_high)
{
  Iflt working_time = t_low;
  Iflt timescale = 1e-10 * length / fL.F_firstDeriv_max(length);
  bool fwdWorking = false;
  
  size_t w = 0;

  while(t_low < t_high)
    {
      //Always try again from the other side
      fwdWorking = !fwdWorking;
    
      if(++w > 1000)
	{
	  std::cerr << "\nWindow shrunk thousands of times\n";
	  
	  return working_time;
	}
    
      working_time = (fwdWorking ? t_low : t_high);
      T tempfL(fL);
      tempfL.stream(working_time);
    
      Iflt deltaT;
      {
	Iflt f0 = tempfL.F_zeroDeriv(),
	  f1 = tempfL.F_firstDeriv(),
	  halff2 = 0.5 * tempfL.F_secondDeriv(),
	  halff2max = 0.5 * tempfL.F_secondDeriv_max(length);
	
	if (f0 > 0) halff2max = -halff2max;
	
	{
	  Iflt boundEnhancer;
	  // Enhance bound, no point continuing if the bounds are out of bounds
	  if (fwdWorking)
	    { if (!quadSolve<ROOT_SMALLEST_POSITIVE>(f0, f1, halff2max, boundEnhancer)) break; }
	  else
	    if (!quadSolve<ROOT_SMALLEST_NEGATIVE>(f0, f1, halff2max, boundEnhancer)) break;
	  
	  (fwdWorking ? t_low : t_high) += boundEnhancer;
	}
	
	if (!quadSolve<ROOT_SMALLEST_POSITIVE>(f0, f1, halff2, deltaT))
	  continue;
      }
      
      if (((working_time + deltaT) > t_high) 
	  || ((working_time + deltaT) < t_low))
	continue;
      
      for(size_t i(1000); i != 0; --i)
	{
	  working_time += deltaT;
	  
	  if((working_time > t_high) || (working_time < t_low))
	    break;
	  
	  tempfL.stream(deltaT);
	  
	  if (!quadSolve<ROOT_SMALLEST_EITHER>(tempfL.F_zeroDeriv(), 
					       tempfL.F_firstDeriv(), 
					       Iflt(0.5 * tempfL.F_secondDeriv()), deltaT))
	    break;
	  
	  if(fabs(deltaT) <  timescale)
	    return working_time + deltaT;
	}
    }

  return HUGE_VAL;
}

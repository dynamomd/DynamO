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

#pragma once

#include "../../../extcode/mathtemplates.hpp"

template<class T>
Iflt quadRootHunter(const T& fL, Iflt length, Iflt& t_low, Iflt& t_high,
		    const Iflt& tolerance)
{
  Iflt working_time = t_low;
  Iflt timescale = tolerance * length / fL.F_firstDeriv_max(length);
  bool fwdWorking = false;

  size_t w = 0;

  while(t_low < t_high)
    {
      //Always try again from the other side
      fwdWorking = !fwdWorking;

      if(++w > 10000)
      {
	std::cerr << "\nThe Frenkel rootfinder is converging too slowly."
	  << "\nt_low = " << t_low << ", t_high = " << t_high;

	if(fabs(t_high - t_low) < timescale)
	{
	  std::cerr << "\nThe gap is small enough to consider the root solved at t_low";
	  return t_low;
	}
	else
	{
	  std::cerr << "\nThe gap is too large and is converging too slowly."
	    << "\n This rootfinding attempt will be aborted and the collision skipped."
	  return HUGE_VAL;
	}
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

	  if (fwdWorking)
	    { if (!quadSolve<ROOT_SMALLEST_POSITIVE>(f0, f1, halff2, deltaT)) continue; }
	  else
	    { if (!quadSolve<ROOT_SMALLEST_NEGATIVE>(f0, f1, halff2, deltaT)) continue; }
	}

      }

      if (((working_time + deltaT) > t_high)
	  || ((working_time + deltaT) < t_low))
	continue;

      // Give it 10,000 iterations before we try shrinking the windows again
      for(size_t i(10000); i != 0; --i)
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

  /* \brief For line line collisions, determines intersections of the infinite lines
  **
  **   Firstly, search for root in main window
  **  - If a root is not found, return failure
  **
  ** If a root is found: bring in an artificial new high boundary just beneath new root
  **  - If this leaves a window, search window for root
  **    - If a root is found, return to top of this section storing only this new root
  **    - If no root is found, drop out of this inner loop
  **  - Check root validity
  **    - If root is valid, this is earliest possible root - roll with it
  **    - If root is invalid, set new concrete t_low just above this found root and go from the top
  */
template<class T>
Iflt frenkelRootSearch(const T& fL, Iflt length, Iflt t_low, Iflt t_high,
		       Iflt tol = 1e-10)
{
  Iflt root = 0.0;

  while(t_high > t_low)
    {
      root = quadRootHunter<T>(fL, length, t_low, t_high, tol);

      if (root == HUGE_VAL) return HUGE_VAL;

      Iflt temp_high = t_high;

      do {
	// Artificial boundary just below root
	T tempfL(fL);
	tempfL.stream(root);

	Iflt Fdoubleprimemax = tempfL.F_secondDeriv_max(length);

	temp_high = root - (fabs(2.0 * tempfL.F_firstDeriv())
			    / Fdoubleprimemax);

	if ((temp_high < t_low) || (Fdoubleprimemax == 0)) break;

	Iflt temp_root = quadRootHunter<T>(fL, length, t_low, temp_high, tol);

	if (temp_root == HUGE_VAL)
	  break;
	else
	    root = temp_root;

      } while(temp_high > t_low);

      // At this point $root contains earliest valid root guess.
      // Check root validity.
      T tempfL(fL);
      tempfL.stream(root);

      if (tempfL.test_root(length))
        return root;
      else
        t_low = root + ((2.0 * fabs(tempfL.F_firstDeriv()))
			/ tempfL.F_secondDeriv_max(length));
    }

  return HUGE_VAL;
}

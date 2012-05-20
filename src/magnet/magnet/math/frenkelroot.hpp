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

#include <magnet/math/quadratic.hpp>
#include <iostream>

namespace magnet {
  namespace math {
    /*! \brief Shooting root finder using quadratic estimation.
  
      \param toleranceLengthScale should be 10^-10 the typical scale of
      the function.
    */
    template<class T>
    std::pair<bool,double> quadRootHunter(const T& fL, double& t_low, double& t_high,
					  const double& toleranceLengthScale)
    {
      double working_time = t_low;
      double timescale = toleranceLengthScale / fL.F_firstDeriv_max();
      bool fwdWorking = false;

      size_t w = 0;

      while(t_low < t_high)
	{
	  //Always try again from the other side
	  fwdWorking = !fwdWorking;

	  if(++w > 1000)
	    {
#ifdef MAGNET_DEBUG
	      std::cerr << "\nThe Frenkel rootfinder is converging too slowly."
			<< "\nt_low = " << t_low << ", t_high = " << t_high;
#endif

	      if(fabs(t_high - t_low) < timescale)
		{
#ifdef MAGNET_DEBUG
		  std::cerr << "\nThe gap is small enough to consider the root solved at t_low";
#endif
		  return std::pair<bool,double>(true, t_low);
		}
	      else
		{
#ifdef MAGNET_DEBUG
		  std::cerr << "\nThe gap is too large and is converging too slowly."
			    << "\n This rootfinding attempt will be aborted and a fake collision returned.";
#endif
		  return std::pair<bool,double>(false, t_low);
		}
	    }

	  working_time = (fwdWorking ? t_low : t_high);
	  T tempfL(fL);
	  tempfL.stream(working_time);

	  double deltaT;
	  {
	    double f0 = tempfL.F_zeroDeriv(),
	      f1 = tempfL.F_firstDeriv(),
	      halff2 = 0.5 * tempfL.F_secondDeriv(),
	      halff2max = 0.5 * tempfL.F_secondDeriv_max();

	    if (f0 > 0) halff2max = -halff2max;

	    {
	      double boundEnhancer;
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
									       double(0.5 * tempfL.F_secondDeriv()), deltaT))
		break;

	      if(fabs(deltaT) <  timescale)
		return std::pair<bool,double>(true, working_time + deltaT);
	    }
	}

      return std::pair<bool,double>(false, HUGE_VAL);
    }

    /*! \brief A root finder that is guarranteed to find the earliest
     root in an interval, for functions with known maximum first and
     second derivatives.
    
     First, search for root in main window
      - If a root is not found, return a failure
    
     If a root is found: start a new search in the window just between this root and the lower bound.
        - If a root is found, restart the search again, in the smaller 
        - If no root is found, drop out of this inner loop
      - Check root validity
        - If root is valid, this is earliest possible root - roll with it
        - If root is invalid, set new concrete t_low just above this found root and go from the top
    
     \param toleranceLengthScale Should be 10^-10 the typical length scale of the system

     If the root finder fails, it will return (false, HUGE_VAL), if it
     fails due to the iteration count going too high, it will return
     (false, t_low), where t_low is a lower bound on any possible root.

     Otherwise it will return (true, t) where t is the location of the
     root.
    */
    template<class T>
    std::pair<bool,double> frenkelRootSearch(const T& fL, double t_low, double t_high,
					     double toleranceLengthScale)
    {
      std::pair<bool,double> root(false,HUGE_VAL);

      while(t_high > t_low)
	{
	  root = quadRootHunter<T>(fL, t_low, t_high, toleranceLengthScale);
	  //M_throw() << "We return " << root.first;
	  //If no root was found, return the lower bound on the root
	  if (root.first == false) return root;

	  //We found a root, now check for earlier roots
	  double temp_high = t_high;
	  do {
	    //Start a search, stream the function to the root
	    T tempfL(fL);
	    tempfL.stream(root.second);
	    //Calculate the offset for the upper bound
	    double Fdoubleprimemax = tempfL.F_secondDeriv_max();
	    temp_high = root.second - (fabs(2.0 * tempfL.F_firstDeriv())
				       / Fdoubleprimemax);

	    //Now check if the upper bound is below the lower bound.
	    //If so, the current root is the earliest.
	    if ((temp_high < t_low) || (Fdoubleprimemax == 0)) break;

	    //Search for a root in the new interval
	    std::pair<bool,double> temp_root = quadRootHunter<T>(fL, t_low, temp_high, toleranceLengthScale);

	    //If there is no root, then the current root is fine
	    if (!temp_root.first)
	      break;
	
	    //Otherwise, use the new root as the earliest one so far
	    root = temp_root;

	  } while(temp_high > t_low);

	  // At this point $root contains earliest valid root guess.
	  // Check root validity.
	  T tempfL(fL);
	  tempfL.stream(root.second);

	  if (tempfL.test_root()) return root;

	  //The root was not valid, set the lower bound to the current root value
	  t_low = root.second + ((2.0 * fabs(tempfL.F_firstDeriv()))
				 / tempfL.F_secondDeriv_max());

	  //Now invalidate the current root
	  root.first = false;
	  root.second = HUGE_VAL;
	}

      return root;
    }
  }
}

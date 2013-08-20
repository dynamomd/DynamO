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
      double timescale = toleranceLengthScale / fL.template max<1>();
      bool fwdWorking = false;

      size_t w = 0;

      while(t_low < t_high)
	{
	  //Always try again from the other side
	  fwdWorking = !fwdWorking;

	  if(++w > 100)
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
	    double f0 = tempfL.template eval<0>(),
	      f1 = tempfL.template eval<1>(),
	      halff2 = 0.5 * tempfL.template eval<2>(),
	      halff2max = 0.5 * tempfL.template max<2>();

	    //This guarantees that the worst-case approximation has
	    //roots on either side of the current time.
	    if (f0 > 0) halff2max = -halff2max;
	    
	    std::pair<double, double> worst_case_roots;
	    try {
	      worst_case_roots = quadraticEquation(halff2max, f1, f0);
	    } catch (NoQuadraticRoots&)
	      {
		M_throw() << "When trying to improve the bounds using worst-case estimates, it was "
		  "found that they could not be improved. This implies there is "
		  "implementation error (zero max 2nd deriv?) in the passed function.";
	      }

	    //Sort the roots
	    if (worst_case_roots.first > worst_case_roots.second) std::swap(worst_case_roots.first, worst_case_roots.second);

#ifdef MAGNET_DEBUG
	    if ((worst_case_roots.first > 0) || (worst_case_roots.second < 0))
	      M_throw() << "The worst case estimates for the root of the function are not either "
		"side of the current location. This implies there is an implementation "
		"error or an untreated numerical edge case.";
#endif
	    //Improve the current boundary
	    if (fwdWorking)
	      t_low += worst_case_roots.second;
	    else
	      t_high += worst_case_roots.first;
	    
	    //Now perform the first step of the shooting
	    std::pair<double, double> estimate_roots;
	    try {
	      estimate_roots = quadraticEquation(halff2, f1, f0);
	    } catch (NoQuadraticRoots&)
	      //If the shooting fails, restart from the other boundary
	      { continue; }
	    
	    //Sort the roots
	    if (estimate_roots.first > estimate_roots.second) std::swap(estimate_roots.first, estimate_roots.second);
	    
	    if (fwdWorking)
	      {
		//Check if there are no positive roots (restart from other boundary if so)
		if (estimate_roots.second < 0) continue;
		//Set deltaT to the smallest positive root.
		if (estimate_roots.first > 0)
		  deltaT = estimate_roots.first;
		else
		  deltaT = estimate_roots.second;
	      }
	    else
	      {
		//Check if there are no negative roots (restart from other boundary if so)
		if (estimate_roots.first > 0) continue;
		//Set deltaT to the smallest negative root.
		if (estimate_roots.second > 0)
		  deltaT = estimate_roots.first;
		else
		  deltaT = estimate_roots.second;
	      }
	  }

	  //Check this first step is still within the other bound
	  if (((working_time + deltaT) > t_high)
	      || ((working_time + deltaT) < t_low))
	    continue;

	  // Give it 100 iterations before we try shrinking the windows again
	  for(size_t i(100); i != 0; --i)
	    {
	      working_time += deltaT;

	      if((working_time > t_high) || (working_time < t_low))
		break;

	      tempfL.stream(deltaT);

	      try {
		std::pair<double, double> estimate_roots = quadraticEquation(0.5 * tempfL.template eval<2>(), tempfL.template eval<1>(), tempfL.template eval<0>());
		if (std::abs(estimate_roots.first) < std::abs(estimate_roots.second))
		  deltaT = estimate_roots.first;
		else
		  deltaT = estimate_roots.second;
	      } catch (NoQuadraticRoots&)
		//If the shooting fails, quit the loop
		{ break; }

	      if (fabs(deltaT) <  timescale)
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
	  //If no root was found, it might have only established a
	  //lower bound, return this lower bound.
	  if (root.first == false) return root;

	  //We found a root, now check for earlier roots in the same interval
	  double temp_high = t_high;
	  do {
	    //Start a search, stream the function to the root
	    T tempfL(fL);
	    tempfL.stream(root.second);
	    //Calculate the offset for the upper bound
	    double Fdoubleprimemax = tempfL.template max<2>();
	    temp_high = root.second - (fabs(2.0 * tempfL.template eval<1>()) / Fdoubleprimemax);

	    //Now check if the upper bound is below the lower bound.
	    //If so, the current root is the earliest.
	    if ((temp_high < t_low) || (Fdoubleprimemax == 0)) break;

	    //Search for a root in the new interval
	    std::pair<bool,double> temp_root = quadRootHunter<T>(fL, t_low, temp_high, toleranceLengthScale);

	    //If there is no root found in the interval, then the current root is fine
	    if ((!temp_root.first) && (temp_root.second == HUGE_VAL)) break;
	
	    //If we have been unable to establish if there is a root
	    //in the interval, we can use the returned time as a lower
	    //bound to come back to later
	    if (!temp_root.first) return temp_root;

	    //Otherwise, the new root is valid, go around and check
	    //for a root in the remaining interval again.
	    root = temp_root;
	  } while(temp_high > t_low);

	  //At this point "root" contains earliest valid root guess.
	  //We now check if this root is acceptable, if not we'll have
	  //to go around again. Fortunately all roots are acceptable
	  //for most algorithms.
	  T tempfL(fL);
	  tempfL.stream(root.second);

	  if (tempfL.test_root()) return root;

	  //The root was not valid, set the lower bound to the current root value
	  t_low = root.second + ((2.0 * fabs(tempfL.template eval<1>())) / tempfL.template max<2>());
	  //and invalidate the current root, resetting the root
	  //finding but with a new lower bound.
	  root.first = false;
	  root.second = HUGE_VAL;
	}

      return root;
    }
  }
}

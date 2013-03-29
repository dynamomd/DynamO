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

#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace {
  const double PI = std::atan(1)*4;
}

namespace dynamo {
  void 
  PotentialLennardJones::outputXML(magnet::xml::XmlStream& XML) const {
    using namespace magnet::xml;
    XML << attr("Type") << "LennardJones"
	<< attr("Sigma") << _sigma
	<< attr("Epsilon") << _epsilon
	<< attr("CutOff") << _cutoff
	<< attr("AttractiveSteps") << _attractiveSteps
      ;
    switch (_U_mode) {
    case MIDPOINT: XML << attr("UMode") << "Midpoint"; break;
    case LEFT: XML << attr("UMode") << "Left"; break;
    case RIGHT: XML << attr("UMode") << "Right"; break;
    case VOLUME: XML << attr("UMode") << "Volume"; break;
    case VIRIAL: 
      XML << attr("UMode") << "Virial"
	  << attr("Temperature") << _kT
	;
      break;
    case MIDVOLUME: XML << attr("UMode") << "MidVolume"; break;
    default:
      M_throw() << "Unknown UMode";
    }

    switch (_R_mode) {
    case DELTAR: XML << attr("RMode") << "DeltaR"; break;
    case DELTAU: XML << attr("RMode") << "DeltaU"; break;
    case DELTAV: XML << attr("RMode") << "DeltaV"; break;
    default:
      M_throw() << "Unknown RMode";
    }
  }

  double 
  PotentialLennardJones::U_uncut(double r) const {
    return 4 * _epsilon * (std::pow(_sigma / r, 12.0) - std::pow(_sigma / r, 6.0));
  }

  double 
  PotentialLennardJones::U(double r) const {
    return U_uncut(r) - U_uncut(_cutoff);
  }

  double 
  PotentialLennardJones::minimum() const {
    return _sigma * std::pow(2.0, 1.0/6.0);
  }

  PotentialLennardJones::PotentialLennardJones(double sigma, double epsilon, double cutoff, UMode umode, RMode rmode, double attractivesteps, double kT):
    _sigma(sigma), _epsilon(epsilon), _cutoff(cutoff), _kT(kT), _attractiveSteps(attractivesteps), _U_mode(umode), _R_mode(rmode)
  {
    _r_cache.push_back(_cutoff);
  }
  
  void 
  PotentialLennardJones::operator<<(const magnet::xml::Node& XML) {
    _r_cache.clear();
    _u_cache.clear();

    _sigma = XML.getAttribute("Sigma").as<double>();
    _epsilon = XML.getAttribute("Epsilon").as<double>();
    _cutoff = XML.getAttribute("CutOff").as<double>();
    
    if (_cutoff <= minimum())
      M_throw() << "The cutoff (" << _cutoff << ") cannot be before the minimum (" 
		<< minimum() << ") in the potential for this Lennard-Jones potential due to the stepping parameters used. Please use a WCA potential instead (if available).";

    _r_cache.push_back(_cutoff);
    
    _attractiveSteps = XML.getAttribute("AttractiveSteps").as<double>();

    const std::string umode_string = XML.getAttribute("UMode").as<std::string>();
    if (!umode_string.compare("Midpoint"))    _U_mode = MIDPOINT;
    else if (!umode_string.compare("Left"))   _U_mode = LEFT;
    else if (!umode_string.compare("Right"))  _U_mode = RIGHT;
    else if (!umode_string.compare("Volume")) _U_mode = VOLUME;
    else if (!umode_string.compare("MidVolume")) _U_mode = MIDVOLUME;
    else if (!umode_string.compare("Virial")) 
      {
	_kT = XML.getAttribute("Temperature").as<double>();
	_U_mode = VIRIAL;
      }
    else
      M_throw() << "Unknown LennardJones UMode (" << umode_string << ") at " << XML.getPath();

    const std::string rmode_string = XML.getAttribute("RMode").as<std::string>();
    if (!rmode_string.compare("DeltaR"))      _R_mode = DELTAR;
    else if (!rmode_string.compare("DeltaU")) _R_mode = DELTAU;
    else if (!rmode_string.compare("DeltaV")) _R_mode = DELTAV;
    else
      M_throw() << "Unknown LennardJones RMode (" << rmode_string << ") at " << XML.getPath();
  }

  std::size_t 
  PotentialLennardJones::steps() const {
    switch (_R_mode) {
    case DELTAR:
      {
	const double deltaR = (_cutoff - minimum()) / _attractiveSteps;
	double steps = _cutoff / deltaR;
	//Rounding down is performed by the conversion to size_t, but
	//we should ensure that any step at r=zero is not
	//included. Just check if the steps variable is a whole integer.
	return size_t(steps) - (size_t(steps) == steps) + 1;
      }
    case DELTAU:
      //In energy stepping there are an infinite number of steps
      return std::numeric_limits<size_t>::max();
    case DELTAV:
      {
	const double deltaV = 4 * PI * (std::pow(_cutoff,3) - std::pow(minimum(),3)) / ( 3 * _attractiveSteps);
	double steps = 4 * PI * std::pow(_cutoff, 3) / (3 * deltaV);
	//Rounding down is performed by the conversion to size_t, but
	//we should ensure that any step at r=zero is not
	//included. Just check if the steps variable is a whole
	//integer.
	return size_t(steps) - (size_t(steps) == steps) + 1;
      }
    default:
      M_throw() << "Unknown RMode";
    }
  }

  long double
  PotentialLennardJones::B2func(const long double r) const {
    return - r * r * (std::exp(-U(r) / _kT) - 1.0);
  }

  void 
  PotentialLennardJones::calculateToStep(const size_t step_id) const {
    const double rmin = minimum();

    //Find the step locations first. we always need one more cached
    //step position than energy, as we need to know the limits of a
    //step to calculate its energy. 
    switch (_R_mode) {
    case DELTAR:
      {
	const double deltaR = (_cutoff - rmin) / _attractiveSteps;

#ifdef DYNAMO_DEBUG
	if (step_id >= steps())
	  M_throw() << "Requested step number " << step_id + 1 << " but there are only " << steps() << " steps in the potential";
#endif

	for (size_t i(_r_cache.size()); i <= step_id; ++i)
	  _r_cache.push_back(_cutoff - i * deltaR);
	
	//Make sure there is one extra step added, and that the zero
	//is added if the end of the stepping is reached
	if (_r_cache.size() == step_id + 1)
	  {
	    if ((step_id == steps() - 1) && (_r_cache.size() == steps()))
	      _r_cache.push_back(0);
	    else
	      _r_cache.push_back(_cutoff - (step_id+1) * deltaR);
	  }
	break;
      }
    case DELTAU:
      {
	const double deltaU = - U(rmin) / _attractiveSteps;
	const size_t minimum_step = size_t(-U(rmin) / deltaU);

	for (size_t i(_r_cache.size()); i <= step_id + 1; ++i)
	  {
	    //Here, we perform a bisection to find the position that
	    //corresponds to an energy. We always perform the
	    //bisection between two limits minR (the r with the lowest
	    //energy) and maxR (the r with the highest energy). We
	    //have to determine these limits and the target U value
	    //depending on if we're before the minimum in the
	    //potential or not.

	    //Assuming we're on a step ID which is before or on the
	    //minimum step. The target energy is decreasing with i
	    //from zero at and the step is bisected by the previous
	    //step and the potential minimum.
	    double target_U = - double(i) * deltaU;
	    double maxR = _r_cache[i-1];
	    double minR = minimum();
	    
	    //Now catch the case if we're after the minimum step
	    //ID. Step upwards in energy with i from the value of the
	    //energy step below the minimum. The last step is the
	    //lower bound for the search but we need to find an upper
	    //limit on r for the bisection. Here just use a expanding
	    //interval which approaches 0.
	    if (i > minimum_step)
	      {
		target_U = (double(i) - 2 * double(minimum_step) - 1) * deltaU;
		minR = std::min(_r_cache[i-1], minimum());
		maxR = minR / 2;
		while (U(maxR) < target_U) maxR /= 2;
	      }


	    //Here we perform the bisection, maxR is above the target
	    //U and minR is below it.
	    for (size_t i(0); i < 1000; ++i)
	      {		  
		double targetR = (maxR + minR) * 0.5;
		double UDiff = U(targetR) - target_U;
		if (UDiff > 0)
		  maxR = targetR;
		else
		  minR = targetR;
		
		if (std::abs(UDiff) <= deltaU * 1e-15) break;
	      }

	    _r_cache.push_back((maxR + minR) * 0.5);
	  }
      }
      break;
    case DELTAV:
      for (size_t i(_r_cache.size()); i <= step_id; ++i)
	_r_cache.push_back(std::pow(std::pow(_r_cache[i-1], 3)- (std::pow(_cutoff,3) - std::pow(minimum(),3)) / _attractiveSteps, 1.0/3.0));

      //Make sure there is one extra step added, and that the zero
      //is added if the end of the stepping is reached
      if (_r_cache.size() == step_id + 1)
	{
	  if ((step_id == steps() - 1) && (_r_cache.size() == steps()))
	    _r_cache.push_back(0);
	  else
	    _r_cache.push_back(std::pow(std::pow(_r_cache[step_id], 3)- (std::pow(_cutoff,3) - std::pow(minimum(),3)) / _attractiveSteps, 1.0/3.0));
	}

      break;
    default:
      M_throw() << "Unknown RMode";
   }

    for (size_t i(_u_cache.size()); i <= step_id; ++i) 
      {
	double newU;
	const double r1 = _r_cache[i+1], r2 = _r_cache[i];
	switch (_U_mode) 
 	  {
	  case MIDPOINT:
	    newU = U((r1 + r2) * 0.5); 
	    break;
	  case LEFT:
	    //Specially treat the case where r=0 is included in the
	    //step.
	    if (r1 == 0)
	      newU = std::numeric_limits<double>::infinity();
	    else
	      newU = U(r1);
	    break;
	  case RIGHT: 
	    newU = U(r2); break;
	  case VOLUME:
	    {
	      const double sigma6 = std::pow(_sigma, 6);
	      const double ri3 = std::pow(r2, 3);
	      const double riplus3 = std::pow(r1, 3);
	      newU = (4 * _epsilon * sigma6 / (ri3 - riplus3)) * (1/ri3 - 1/riplus3 - (sigma6/3.0) * (1/(ri3*ri3*ri3) - 1/(riplus3*riplus3*riplus3))) - U_uncut(_cutoff);

	      //Specially treat the case where r=0 is included in the
	      //step. The singularity dominates and the step energy is
	      //infinite.
	      if (r1 == 0) newU = std::numeric_limits<double>::infinity();
	    }
	    break;
	  case VIRIAL:
	    {
	      //Numerically integrate for the B2 in the region [r1,r2]
	      //using simpsons rule

	      //Must be even
	      const size_t iterations = 100000;
	      const long double h= (r2 - r1) / (long double)(iterations);
	      long double sum(0);

	      //Specially treat the case where r=0 is included in the
	      //step. It doesn't contribute to the virial provided the
	      //temperature is finite, but the calculation has a
	      //divide by zero in it so it must be avoided.
	      if (r1 != 0) sum += B2func(r1);

	      for (size_t i(1); i < iterations; ++i)
		if (i % 2)
		  sum += 4 * B2func(r1 + i * h);
		else
		  sum += 2 * B2func(r1 + i * h);
	      
	      sum += B2func(r2);

	      long double B2 = h * sum / 3;

	      long double log_arg = 1 - 3 * B2/ (r2 * r2 * r2 - r1 * r1 * r1);
	      newU = - _kT * std::log(log_arg);

	      //Sometimes, precision errors in the log cause a
	      //negative log_arg, if this is the case the potential
	      //energy is effectively infinity.
	      if (log_arg <= 0) newU = std::numeric_limits<double>::infinity();
	    }
	    break;
	  case MIDVOLUME:
	    {
	      newU = U(std::pow((r1*r1*r1 + r2*r2*r2)/2, 1.0/3.0)); 
	      break;
	    }
	  default: M_throw() << "Unknown UMode";
	  }
	_u_cache.push_back(newU);
      }
  }
}
